function result = simulateEcholocation(bandwidth, kr, target_distance, initial_velocity, initial_call_duration, motile, makeAudio)
  
% simulateEcholocation
% Simulates an FM bat's echolocation behaviour while approaching a target. 
% The model includes call timing, target motion, amplitude adaptation, and echo simulation.
% USAGE:
%   result = simulateEcholocation(bandwidth, kr, target_distance, ...
%                                 initial_velocity, initial_call_duration, ...
%                                 motile, makeAudio)
%
% INPUTS:
%   bandwidth             - [Hz] 1x2 vector [minFreq, maxFreq] defining the frequency range of the call
%   kr                    - [unitless] call rate responsivity coefficient (scales echo delay to next call delay)
%   target_distance       - [m] initial distance to the target
%   initial_velocity      - [m/s] velocity of approach (positive toward target)
%   initial_call_duration - [s] duration of calls at the beginning
%   motile                - [logical] whether the bat changes its velocity due to target motion (adds noise)
%   makeAudio             - [logical] if true, generates synthetic audio sequences for calls and echoes
%
% OUTPUT:
%   result                - struct containing:
%       .delta_s               - stepwise distance changes between calls
%       .delta_t               - inter-pulse intervals (IPI)
%       .call_rate_responsivity - inverse of IPI (Hz)
%       .call_levels           - calibrated call amplitudes (Pa)
%       .call_points           - call timestamps (sample indices)
%       .decrement_point       - index at which duration starts decreasing
%       .detection_point       - index at which echo detection drops below threshold
%       .echo_points           - echo return timestamps (sample indices)
%       .echo_time             - time of each echo arrival (s)
%       .initial_velocity      - input parameter, for reference
%       .new_amplitude         - updated amplitude per call
%       .new_call_duration     - updated call duration array
%       .initial_target_distance - input parameter, for reference
%       .new_target_distance   - updated target distance after each call
%       .Ta                    - two-way travel time per call (s)
%       .Tb                    - inter-pulse interval at max responsivity
%       .Tb_prime              - rate of change of IPI near Tb
%       .Tb_idx                - index corresponding to Tb
%       .time_to_contact       - estimated remaining time to reach target
%       .TS                    - target strength per call
%       .Kr                    - input parameter, for reference
%       .max_call_rate         - max physiologically plausible call rate (Hz)
%       .MaxCallRatePoint      - index when IPI exceeds max rate
%       .OverlapPoint          - index where call-echo overlap starts
%       .audio                 - 2-column audio waveform [call, echo] if makeAudio is true; empty otherwise
%
% NOTES:
% - This model assumes free-field spherical propagation and adjusts signal level based on geometric and air attenuation.
% - Target motion can include noise for more realism (when `motile` is true).
% - Audio synthesis uses `generateVirtualBatCall` and `airAttenuationFilter`.
% - Echoes are created by convolving the call with a distance-dependent air attenuation filter.
%
% EXAMPLE:
%   bw = [45000 90000]; % 45–90 kHz
%   kr = 0.5;
%   d = 1.0;            % 1 metre
%   v = 2.0;            % 2 m/s approach
%   duration = 0.003;   % 3 ms
%   motile = true;
%   makeAudio = true;
%   result = simulateEcholocation(bw, kr, d, v, duration, motile, makeAudio);
    fs = 192e3;
    c = 343;
    target_diameter = freq2wavelen(min(bandwidth), c); % m
    source_level = 94; % dB
    source_level_minimum = 70; % dB
    SNR = 20; % dB
    detection_level = 20; % dB
    max_call_rate = 200; % calls/second
    temperature = 25; % °C
    humidity = 50; % %

    alfa = atmatt(temperature, humidity, bandwidth);
    A = alfa(bandwidth == max(bandwidth));

    new_target_distance(1) = target_distance;
    Ta(1) = 2 * target_distance / c;
    delta_t(1) = 0;
    echo_time(1) = delta_t(1) + Ta(1);
   

    call_number = 1;
    done = false;
    while ~done
        Ta(call_number + 1) = 2 * new_target_distance(call_number) / c;
        delta_t(call_number + 1) = kr * Ta(call_number + 1);
        echo_time(call_number + 1) = max(cumsum(delta_t(1:call_number + 1))) + Ta(call_number + 1);
        if motile
            delta_s(call_number) = delta_t(call_number + 1) * initial_velocity + randn() / 10;
        else
            delta_s(call_number) = delta_t(call_number + 1) * initial_velocity; % + randn() / 10;
        end
        new_target_distance(call_number + 1) = new_target_distance(call_number) - delta_s(call_number);

        TS(call_number) = calculateTargetStrength(min(bandwidth), target_diameter, new_target_distance(call_number + 1), A);
        TL(call_number) = 2 * calculateGeometricAttenuation(new_target_distance(call_number + 1));
        RL(call_number) = source_level - abs(TL(call_number)) + TS(call_number);

        if RL(call_number) < detection_level
            detection_point = call_number + 1;
        end
        if Ta(call_number) > kr * initial_call_duration
            decrement_point = call_number + 1;
        end
        if new_target_distance(call_number) < 0.05
            done = true;
        end
        call_number = call_number + 1;
    end

    call_rate_responsivity = abs(1 ./ diff(delta_t));
   
    new_call_duration = initial_call_duration * ones(call_number, 1)';
    new_call_duration(decrement_point:end) = Ta(decrement_point:end) ./ kr;
    new_call_duration(new_call_duration < 0.0005) = 0.0005;

    new_amplitude = source_level * ones(call_number - 1, 1);
    new_amplitude(detection_point:end) = source_level - abs((detection_level - RL(detection_point:end)));
    new_amplitude(new_amplitude < source_level_minimum) = linspace(source_level_minimum, source_level_minimum - 6, length(new_amplitude(new_amplitude < source_level_minimum)));
    call_levels = 20e-6 * 10.^(new_amplitude ./ 20);

    seq = randn(fs * round(max(cumsum(delta_t))), 1) ./ 10e5;
    call_points = round(fs .* cumsum(delta_t));
    echo_seq = randn(round(fs * max(echo_time)), 1) ./ 10e5;
    echo_points = round(fs .* echo_time);

    % calculate call rate
    timestamps = (call_points)./fs;
    call_rate_profile = 1./diff(timestamps);

    overlap_point = find(call_points + round(fs * new_call_duration) > echo_points, 1, 'first');
    overlap_distance = new_target_distance(overlap_point);
    overlap_call_point = call_points(overlap_point);
    max_call_rate_point = find(call_rate_profile > max_call_rate, 1, 'first');
    [~, Tb_idx] = min(abs(call_rate_responsivity - max_call_rate));
    Tb = delta_t(Tb_idx);
    ddt = diff(delta_t);
    Tb_prime = ddt(Tb_idx); % this should be it!
    time_to_contact = new_target_distance/initial_velocity; % remaining reaction time window 

    if makeAudio
        for i = 1:call_number - 1
            call = generateVirtualBatCall(min(bandwidth), max(bandwidth), new_call_duration(i), fs, 0);
            call = call .* call_levels(i);
            [ir, ~] = airAttenuationFilter(new_target_distance(i), fs, 64, 2);
            attenuation = 2 * calculateGeometricAttenuation(new_target_distance(i));
            att_factor = call_levels(i) - 20e-6 * 10^(abs(attenuation) / 20);
            echo = att_factor .* conv(call, ir, "same");
            seq(call_points(i)+1:call_points(i)+length(call)) = call;
            echo_seq(echo_points(i)+1:echo_points(i)+length(echo)) = echo;
        end

        audio_out = [seq(1:min(end, length(echo_seq))), echo_seq(1:min(end, length(seq)))];
    else
        audio_out = [];
    end
    result = struct(...
        'delta_s', delta_s, ...
        'delta_t', delta_t, ...
        'call_rate_responsivity', call_rate_responsivity, ...
        'call_levels', call_levels, ...
        'call_points', call_points, ...
        'decrement_point', decrement_point, ...
        'detection_point', detection_point, ...
        'echo_points', echo_points, ...
        'echo_time', echo_time, ...
        'initial_velocity', initial_velocity,...
        'new_amplitude', new_amplitude, ...
        'new_call_duration', new_call_duration, ...
        'initial_target_distance', target_distance, ...
        'new_target_distance', new_target_distance, ...
        'Ta', Ta, ...
        'Tb', Tb, ... % IPI at at Responsivity max
        'Tb_prime', Tb_prime, ... % diff in IPI at Responsivity max
        'Tb_idx', Tb_idx, ...
        'time_to_contact', time_to_contact, ...
        'TS', TS, ...
        'Kr', kr,...
        'max_call_rate', max_call_rate, ...
        'MaxCallRatePoint', max_call_rate_point,...
        'OverlapPoint', overlap_point,...
        'audio', audio_out ...
    );
end