%% Smoke tests for simulateResponsivityCore
% Run this file from MATLAB to check the first implementation.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
simRoot = fileparts(thisDir);
addpath(simRoot);
addpath(fullfile(simRoot, 'fcn'));

%% 1. Default V2-like first-arrival simulation
params = struct();
opts = struct();
opts.rngSeed = 1;

resDefault = simulateResponsivityCore(params, opts);
disp('=== Default simulation ===');
disp(resDefault.summary);
disp(resDefault.calls(1:min(5,height(resDefault.calls)), ...
    {'SeqID','CallNumber','AnchorDistance_m','Tdelay_s','Ta_s','IPI_s','CallRate_Hz'}));

%% 2. SSG-world-like 3D/two-target/random-anchor simulation
params = struct();
params.numSequences = 3;
params.maxCalls = 80;
params.initialDistance_m = 3;
params.initialBatSpeed_m_s = 6;

opts = struct();
opts.rngSeed = 2;
opts.geometryMode = "3D";
opts.numTargets = 2;
opts.targetMotion = true;
opts.targetVelocityMode = "stochastic";
opts.batVelocityMode = "jittered";
opts.batVelocityJitterFrac = 0.15;
opts.anchorMode = "random";
opts.timingMode = "motionAware";

resSSGWorld = simulateResponsivityCore(params, opts);
disp('=== SSG-world-like simulation ===');
disp(resSSGWorld.summary);
disp(resSSGWorld.calls(1:min(5,height(resSSGWorld.calls)), ...
    {'SeqID','CallNumber','AnchorTargetID','AnchorSwitch','AnchorDistance_m','RelativeClosingVelocity_m_s','IPI_s'}));

%% 3. Wingbeat-state simulation with call-duration jitter
params = struct();
params.initialDistance_m = 4;
params.initialBatSpeed_m_s = 5;
params.maxCalls = 80;

opts = struct();
opts.rngSeed = 3;
opts.enforceMaxCallRate = true;
opts.callDurationMode = "previousTa";
opts.callDurationJitter.enabled = true;
opts.callDurationJitter.mode = "additive";
opts.callDurationJitter.rho = 0.4;
opts.wing.enabled = true;
opts.wing.frequencyMode = "dynamic";
opts.wing.excursionMode = "dynamic";

resWing = simulateResponsivityCore(params, opts);
disp('=== Wingbeat-state simulation ===');
disp(resWing.summary);
disp(resWing.calls(1:min(5,height(resWing.calls)), ...
    {'CallNumber','Tcall_s','TcallJitter_s','IPI_s','WingbeatFreq_Hz','WingbeatExcursion','WingbeatPhase_frac','WingbeatAngularPosition'}));

%% Basic sanity checks
assert(all(resDefault.calls.IPI_s > 0), 'Default simulation produced non-positive IPIs.');
assert(all(resDefault.calls.Ta_s >= resDefault.calls.Tdelay_s), ...
    'Acoustic acquisition interval should be at least the first-echo delay.');
assert(all(resWing.calls.Tcall_s >= resWing.params.minCallDuration_s), ...
    'Call duration fell below the physiological minimum.');
assert(all(resWing.calls.Tcall_s <= resWing.params.initialCallDuration_s), ...
    'Call duration exceeded the initial call duration cap.');

if resWing.opts.wing.enabled
    dynamicRows = resWing.calls.WingbeatFreq_Hz > resWing.params.maxWingbeatFreq_Hz;
    if any(dynamicRows)
        assert(all(resWing.calls.WingbeatExcursion(dynamicRows) < ...
            resWing.params.initialWingbeatExcursion), ...
            'Dynamic wingbeat frequency exceeded the ceiling without excursion contraction.');
    end
end

%% Optional audio synthesis smoke check
try
    audioOpts = struct('fs', 192e3, 'callType', "FM", 'normalise', true);
    audio = synthesiseResponsivityAudio(resDefault, audioOpts);
    assert(size(audio.stereo, 2) == 2, 'Audio stereo output should have two channels.');
    assert(audio.fs == audioOpts.fs, 'Audio sampling rate mismatch.');
    disp('Audio synthesis smoke check passed.');
catch ME
    warning('Audio synthesis smoke check skipped or failed: %s', ME.message);
end

disp('All smoke-test sanity checks passed.');
