function result = simulateEcholocationCF(frequency, kr, target_distance, initial_velocity, initial_call_duration)
    fs = 192e3;
    c = 343;
    source_level = 94; % dB
    source_level_minimum = 70; % dB
    SNR = 20; % dB
    detection_level = 20; % dB
    max_call_rate = 100; % calls/second
    temperature = 25; % °C
    humidity = 50; % %

    % Echo attenuation model
    alfa = atmatt(temperature, humidity, frequency);
    A = alfa;

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
        delta_s(call_number) = delta_t(call_number + 1) * initial_velocity;
        new_target_distance(call_number + 1) = new_target_distance(call_number) - delta_s(call_number);

        TS(call_number) = calculateTargetStrength(frequency, ...
                            freq2wavelen(frequency, c), new_target_distance(call_number + 1), A);
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
    new_amplitude(new_amplitude < source_level_minimum) = ...
        linspace(source_level_minimum, source_level_minimum - 6, length(new_amplitude(new_amplitude < source_level_minimum)));
    call_levels = 20e-6 * 10.^(new_amplitude ./ 20);

    seq = randn(fs * round(max(cumsum(delta_t))), 1) ./ 10e5;
    call_points = round(fs .* cumsum(delta_t));
    echo_seq = randn(round(fs * max(echo_time)), 1) ./ 10e5;
    echo_points = round(fs .* echo_time);

    timestamps = (call_points)./fs;
    call_rate_profile = 1 ./ diff(timestamps);

    overlap_point = find(call_points + round(fs * new_call_duration) > echo_points, 1, 'first');
    overlap_distance = new_target_distance(overlap_point);
    overlap_call_point = call_points(overlap_point);
    [~, Tb_idx] = min(abs(call_rate_responsivity - max_call_rate));
    Tb = delta_t(Tb_idx);

    for i = 1:call_number - 1
        % Emit CF call (no Doppler compression for emission)
        call = generateCFBatCall(frequency, new_call_duration(i) * 1000, fs, 0, 0);
        call = call .* call_levels(i);

        % Echo with Doppler compression
        echo_raw = generateCFBatCall(frequency, new_call_duration(i) * 1000, fs, 0, initial_velocity);
        echo_raw = echo_raw .* call_levels(i);

        [ir, ~] = airAttenuationFilter(new_target_distance(i), fs, 64, 2);
        attenuation = 2 * calculateGeometricAttenuation(new_target_distance(i));
        att_factor = call_levels(i) - 20e-6 * 10^(abs(attenuation) / 20);
        echo = att_factor .* conv(echo_raw, ir, "same");

        seq(call_points(i)+1:call_points(i)+length(call)) = call;
        echo_seq(echo_points(i)+1:echo_points(i)+length(echo)) = echo;
    end

    audio_out = [seq(1:min(end, length(echo_seq))), echo_seq(1:min(end, length(seq)))];

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
        'Tb', Tb, ...
        'Tb_idx', Tb_idx, ...
        'TS', TS, ...
        'Kr', kr,...
        'max_call_rate', max_call_rate, ...
        'MaxCallRatePoint', find(call_rate_profile > max_call_rate, 1, 'first'), ...
        'OverlapPoint', overlap_point, ...
        'audio', audio_out ...
    );
end