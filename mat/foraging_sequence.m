%% Foraging Call Sequence Simulation - Vespertillio numeralis
% -------------------------------------------------------------------------
% This script simulates a biologically-inspired echolocation behaviour 
% during prey pursuit using a closed-loop model based on Responsivity Theory.
%
% Species: Vespertillio numeralis (simulated)
% Model: CRPC (Closed-loop Responsivity-based Pursuit Control)
% Author: Ravi Umadi
% Date: [Insert Date]
%
% The script:
% - Generates a synthetic sequence of bat calls and echoes
% - Updates call parameters based on echo timing and distance
% - Simulates a motile prey target
% - Calculates adjusted call rate, amplitude, and duration
% - Includes attenuation, target strength, and responsivity dynamics
% - Plots waveforms and behaviourally relevant transitions
% - Outputs stereo audio (.wav) of call-echo sequence
%
% -------------------------------------------------------------------------
% Key Parameters:
%   fs                    : Sampling frequency (Hz)
%   bandwidth             : Call bandwidth (start to end frequency in Hz)
%   c                     : Speed of sound (m/s)
%   kr                    : Responsivity coefficient (dimensionless)
%   target_distance       : Initial target distance (m)
%   target_diameter       : Derived from frequency and wavelength (m)
%   initial_velocity      : Bat velocity (m/s)
%   initial_call_duration : Starting duration of each call (s)
%   source_level          : Initial source level of calls (dB)
%   SNR                   : Desired signal-to-noise ratio (dB)
%   max_call_rate         : Cap on call rate (calls/sec)
%   temperature, humidity : For atmospheric attenuation calculation
%
% Outputs:
%   - Time series plot of calls and echoes
%   - Behavioural transitions: detection point, duration decrement, buzz onset
%   - Velocity profile during pursuit
%   - Stereo WAV file: 'seq_3ms.wav'
%
% Dependencies:
%   - freq2wavelen() (phased array toolbox), atmatt(), calculateTargetStrength()
%   - calculateGeometricAttenuation(), airAttenuationFilter()
%   - generateVirtualBatCall(), formatLatex()
%
% -------------------------------------------------------------------------
% Usage Notes:
%   - Designed for ultrasonic simulation (≥192kHz)
%   - Simulates prey movement using added noise in target distance
%   - echo_time and delta_t dynamically control loop timing
%   - This model allows exploration of CRPC-derived behavioural patterns
%     like call rate saturation, echo-call overlap, and buzz phase onset
%
% Example Analysis:
%   - Call rate calculated using diff(delta_t)
%   - Overlap point identified where call ends before echo returns
%   - Tb and Tb' values derived as critical indicators of transition
%
% Citation:
%   This simulation is part of a theoretical exploration of bat biosonar
%   dynamics, building on the concept of sensory–motor responsivity.
%
% -------------------------------------------------------------------------
clear, close
% Parameters
fs = 192e3;
bandwidth = 65e3:-100:35e3;
c = 343;
kr = 5;
target_distance = 20; %m
target_diameter = freq2wavelen(min(bandwidth), c); %m
initial_velocity = 8; %m/s
initial_call_duration = 0.005; % s
source_level = 94; % dB
source_level_minimum = 70; %dB
SNR = 20; % dB
detection_level = 20; % dB
max_call_rate = 200; % calls/second
% buzz_rate = 200;
% atmospheric absorption coefficient for the bandwidth
temperature = 25; % Deg C
humidity = 50; %
alfa = atmatt(temperature,humidity,bandwidth);
A = alfa(bandwidth == max(bandwidth));

%
new_target_distance(1) = target_distance;
Ta(1) = 2*target_distance/c;

delta_t(1) = 0;
echo_time(1) = delta_t(1) + Ta(1);

call_number = 1;
done = false;
while ~done
    Ta(call_number+1) = 2*new_target_distance(call_number)/c;
    delta_t(call_number+1) = (kr*Ta(call_number+1));
    echo_time(call_number+1) = max(cumsum(delta_t(1:call_number+1))) + Ta(call_number+1);
    delta_s(call_number) = delta_t(call_number+1)*initial_velocity + randn()/10; %% add 'randn()/10' for simulating a motile prey. 
    new_target_distance(call_number+1) = new_target_distance(call_number) - delta_s(call_number); % rand to simulate motile prey

    % get the RL
    TS(call_number) = calculateTargetStrength(min(bandwidth), target_diameter, new_target_distance(call_number+1), A);
    TL(call_number) = 2*calculateGeometricAttenuation(new_target_distance(call_number+1));
    RL(call_number) = source_level - abs(TL(call_number)) + TS(call_number);

    % mark detection point, start level adjustment from here
    if RL(call_number) < detection_level
        detection_point = call_number+1;
    end
    % mark call duration decrement point
    if Ta(call_number) > kr*initial_call_duration
        decrement_point = call_number+1;
    end

    if new_target_distance(call_number) < 0.05
        done = true;
    end
    call_number = call_number+1;
    disp(['call ', num2str(call_number)])

    % pause,
end

call_rate_responsivity = abs(1./diff(delta_t)); 



% calculate adjusted call duration
new_call_duration = initial_call_duration* ones(call_number,1)';
new_call_duration(decrement_point:end) = Ta(decrement_point:end)./kr;
new_call_duration(new_call_duration<0.0005) = 0.0005;

% calculate adjusted levels
new_amplitude = source_level * ones(call_number-1,1);
new_amplitude(detection_point:end) = source_level- abs((detection_level - RL(detection_point:end)));
new_amplitude(new_amplitude<source_level_minimum) = linspace(source_level_minimum, source_level_minimum-6, length(new_amplitude(new_amplitude<source_level_minimum))) ;
call_levels = 20e-6* 10.^(new_amplitude./20);

% create the call sequence
seq = randn(fs*round(max(cumsum(delta_t))), 1)./10e5;
call_points = fs.*cumsum(delta_t);
echo_seq = randn(round(fs*max(echo_time)), 1)./10e5;
echo_points = fs.*(echo_time);

% calculate call rate
timestamps = (call_points)./fs;
call_rate_profile = 1./diff(timestamps);

% find overlap points
overlap_point = find(call_points + (fs*new_call_duration) > echo_points, 1, 'first');
overlap_distance = new_target_distance(overlap_point);
overlap_call_point = call_points(overlap_point);
max_call_rate_point = find(call_rate_profile < max_call_rate, 1, 'last');

% Tb_minimum = delta_t(max_call_rate_point) - 0.0005;
[~, Tb_idx] = min(abs(call_rate_responsivity - max_call_rate));
Tb = delta_t(Tb_idx);
ddt = diff(delta_t);
Tb_prime = ddt(Tb_idx); % this should be it!
time_to_contact = new_target_distance/initial_velocity; % remaining reaction time window 
buzz_ready_point = call_points(Tb_idx);
freq_drop_window = 1./(1:call_number);

for i = 1: call_number-1
    
    if i > max_call_rate_point
        call = generateVirtualBatCall(min(bandwidth)-(min(bandwidth)*10*freq_drop_window(i)), max(bandwidth)-max(bandwidth)*10*freq_drop_window(i), new_call_duration(i), fs, 0);
    else
        call = generateVirtualBatCall(min(bandwidth), max(bandwidth), new_call_duration(i), fs, 0);
    end

    % adjust the level
    call = call.*call_levels(i);
    [ir, geom_atten] = airAttenuationFilter(new_target_distance(i), fs, 64, 2);
    attenuation = 2*calculateGeometricAttenuation(new_target_distance(i));
    att_factor = call_levels(i) - 20e-6*10^(abs(attenuation)/20);
    % ir = ir./max(ir);
    echo = att_factor.* conv(call,ir, "same");
    seq(call_points(i)+1:call_points(i)+length(call)) = call;
    echo_seq(echo_points(i)+1:echo_points(i)+length(echo)) = echo;
end


% plot
colors = [
    0.0500 0.05 0.050;      % black
    1.0000 0.5490 0.0000;   % orange
    0.0000 0.4470 0.7410;   % blue
    0.0000 0.6000 0.2000;   % green
    0.8500 0.1250 0.0980    % red
];
linwidth = 2.5;
figure, plot((1:length(seq))./fs, seq, 'k')
hold on
plot((1:length(echo_seq))./fs, echo_seq, 'm')
xline(call_points(detection_point)/fs, 'LineWidth',linwidth, 'Color', colors(1,:), 'LineStyle','-.');
xline(call_points(decrement_point)/fs, 'LineWidth',linwidth, 'Color', colors(2,:), 'LineStyle','-.');
xline(buzz_ready_point/fs, 'LineWidth',linwidth, 'Color', colors(3,:), 'LineStyle','-.');
xline(call_points(max_call_rate_point)/fs, 'LineWidth',linwidth, 'Color', colors(4,:), 'LineStyle','-.');
xline(overlap_call_point/fs, 'LineWidth',linwidth, 'Color', colors(5,:), 'LineStyle','-.')

ylim([-0.1 0.1])
% xlim([0.8, 1.8])
ymax = max(abs(ylim));        % Get the max absolute value for symmetric range
nTicks = 3;                   % Set how many ticks you want (odd number preferred for zero in middle)
yticks(linspace(-ymax, ymax, nTicks));
yt= yticks;
xt = xticks;
% xticklabels(target_distance-(xt*initial_velocity))
xticklabels(sort(xt, 'descend'))
% xlim([4 6.8])
db_scale = real(20.*log10(yt./2e-05));
db_scale(db_scale <0) = 0;
yticklabels(round(db_scale))

% axisbold
% grid minor
% grid on
%
title(['Velocity ' num2str(initial_velocity) 'm/s -- Number of Calls ' num2str(call_number) ' - Time ' num2str(round(max(timestamps),2)) ' s ' '$\& C_{r} ~(@ T_{b^*} = ~$' num2str(abs(round(Tb_prime*1000, 2))) '~ms) - ' num2str(round(call_rate_profile(Tb_idx))) ' Hz'], 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight','bold');
% xlabel('Distance to Target, m', 'Interpreter','latex', 'FontSize', 16, 'FontWeight','bold')
xlabel('Time, s', 'Interpreter','latex', 'FontSize', 16, 'FontWeight','bold')
ylabel('dB re. 20$\mu$Pa', 'Interpreter','latex', 'FontSize', 16, 'FontWeight','bold')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
set(gcf,'units','inch','position',[0,0,18,2]); 
figName = ['foraging_sequence_motile_' num2str(initial_velocity) 'ms'];
figpath = fullfile('', figName); % set your path
formatLatex(gca)
% exportgraphics(gcf, [figpath '.png'], 'Resolution', 300, 'Append', false)
%% optionally plot velocity profile
figure, plot(call_points(2:end)./fs, abs(delta_s./delta_t(2:end)), '-o', 'MarkerFaceColor', [0.2 0.7 0.2]);
    ylabel('Velocity (m/s)', 'Interpreter','latex');
    xlabel('Time (s)', 'Interpreter','latex');
    title('Estimated or Measured Velocity', 'Interpreter','latex', 'FontSize', 14);
    axis tight;
    formatLatex(gca)
    set(gcf,'units','inch','position',[0,0,18,2]); 
%% overlap plots
figure, plot(delta_t(2:end), call_rate_profile)
ylim([0 500])
xlim([0 0.5])
xline(delta_t(overlap_point),'LineWidth',0.75)
xline(delta_t(max_call_rate_point), 'LineWidth',0.75, 'Color', 'r')


%% write to wav file
% format
seq_length_diff = length(seq) - length(echo_seq);
if seq_length_diff > 1
echo_seq = [echo_seq ; randn(seq_length_diff,1)./10e5];
elseif seq_length_diff < 0
    seq = [seq ; randn(seq_length_diff,1)./10e5];
end

out = [seq echo_seq];
out = [randn(fs,2)./10e5; out; randn(fs,2)./10e5];
out = out .*0.95;
audiowrite('seq_3ms.wav', out, fs)