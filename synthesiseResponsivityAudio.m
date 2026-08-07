function audio = synthesiseResponsivityAudio(result, audioOpts)
% synthesiseResponsivityAudio
% Generate call/echo audio from a simulateResponsivityCore result.
%
% This function is intentionally separate from the core simulator. It reads
% the call table and synthesises waveform channels for inspection or demos.
%
% Usage:
%   audio = synthesiseResponsivityAudio(result)
%   audio = synthesiseResponsivityAudio(result, audioOpts)
%
% audioOpts fields:
%   fs              default 192e3
%   callType        "fromTable", "FM", or "CF"
%   normalise       true/false
%   echoGain        scalar echo gain, default 0.25
%   distanceScaling true/false, default true
%   tailPercent     zero padding passed to call generators, default 0

    if nargin < 2 || isempty(audioOpts)
        audioOpts = struct();
    end
    audioOpts = setDefault(audioOpts, 'fs', 192e3);
    audioOpts = setDefault(audioOpts, 'callType', "fromTable");
    audioOpts = setDefault(audioOpts, 'normalise', true);
    audioOpts = setDefault(audioOpts, 'echoGain', 0.25);
    audioOpts = setDefault(audioOpts, 'distanceScaling', true);
    audioOpts = setDefault(audioOpts, 'tailPercent', 0);

    thisDir = fileparts(mfilename('fullpath'));
    fcnDir = fullfile(thisDir, 'fcn');
    if exist(fcnDir, 'dir')
        addpath(fcnDir);
    end

    T = result.calls;
    fs = audioOpts.fs;

    if isempty(T)
        audio = struct('fs', fs, 'call', [], 'echo', [], 'stereo', [], 'opts', audioOpts);
        return
    end

    endTime = max([T.NextCallOnsetTime_s; T.EchoAnchorTime_s]) + max(T.Tcall_s);
    nSamples = max(1, ceil(endTime * fs) + 1);
    callChannel = zeros(nSamples, 1);
    echoChannel = zeros(nSamples, 1);

    for i = 1:height(T)
        callWave = makeCallWave(T(i, :), fs, audioOpts);
        if isempty(callWave)
            continue
        end

        callStart = max(1, round(T.CallOnsetTime_s(i) * fs) + 1);
        callChannel = addWave(callChannel, callWave, callStart);

        echoWave = callWave;
        echoGain = audioOpts.echoGain;
        if audioOpts.distanceScaling
            echoGain = echoGain / max(T.AnchorDistance_m(i)^2, 0.01);
        end
        echoWave = echoGain * echoWave;

        echoStart = max(1, round(T.EchoStartTime_s(i) * fs) + 1);
        echoChannel = addWave(echoChannel, echoWave, echoStart);
    end

    if audioOpts.normalise
        peak = max(abs([callChannel; echoChannel]));
        if peak > 0
            callChannel = callChannel ./ peak;
            echoChannel = echoChannel ./ peak;
        end
    end

    audio = struct();
    audio.fs = fs;
    audio.call = callChannel;
    audio.echo = echoChannel;
    audio.stereo = [callChannel echoChannel];
    audio.opts = audioOpts;
end

function wave = makeCallWave(row, fs, audioOpts)
    callType = string(audioOpts.callType);
    if callType == "fromTable"
        callType = string(row.CallType);
    end

    switch upper(char(callType))
        case 'FM'
            f0 = row.StartFreq_Hz;
            f1 = row.EndFreq_Hz;
            wave = generateVirtualBatCall(f0, f1, row.Tcall_s, fs, audioOpts.tailPercent);
        case 'CF'
            wave = generateCFBatCall(row.CFreq_Hz, 1000 * row.Tcall_s, fs, ...
                audioOpts.tailPercent, row.RelativeClosingVelocityForTiming_m_s);
        otherwise
            error('Unknown call type for synthesis: %s', callType);
    end

    wave = wave(:);
end

function y = addWave(y, wave, startIdx)
    stopIdx = min(numel(y), startIdx + numel(wave) - 1);
    if startIdx > stopIdx
        return
    end
    n = stopIdx - startIdx + 1;
    y(startIdx:stopIdx) = y(startIdx:stopIdx) + wave(1:n);
end

function s = setDefault(s, fieldName, value)
    if ~isfield(s, fieldName) || isempty(s.(fieldName))
        s.(fieldName) = value;
    end
end
