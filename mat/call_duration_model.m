% call duration decrement model
% model: t_call = (T_a + t_b') - \Delta t

c = 343;
target_distance = linspace(20,0.02, 1000)'; % m
subject_velocity = 5; % m/s, since c >> v, it does not matter
call_durations = (1:1:10)*10e-4; % s
kr = (1:10)';
decrement_start_distance = (call_durations.*c)./2;

% Loop over target_distance
for i = 1:length(target_distance)
    % Calculate deltaS, signal_flight_time, and new_target_distance
    [deltaS, signal_flight_time, new_target_distance] = ...
        emissionReceptionSpaceTimeChange(subject_velocity, ...
        target_distance(i));

    % Store results in arrays
    deltaS_results(i) = deltaS;
    signal_flight_time_results(i) = signal_flight_time;
    new_target_distance_results(i) = new_target_distance;
end

% Calculate time-of-flight (TOF) for each distance
Ta = (2 * new_target_distance_results) / c;

% threshhold ta
for i = 1: length(call_durations)
    sd = find(Ta < call_durations(i));
    for j = 1:length(kr)
        new_call_duration = kr(j)*Ta(sd);
        decrement_distance(i,j) = target_distance(find(Ta < max(new_call_duration), 1, "first"));
    end
end

%% plot
% color each line from green to red
color_range = linspace(0, 1, length(call_durations));
fig = figure;
for i = 1: length(call_durations)
    plot(decrement_distance(:,i), call_durations*10e2, 'Color', [color_range(i), 1-color_range(i), 0], 'LineWidth',0.75)
    hold on
    drawnow

end
title('Distance Threshold as function of $t_{call}$ and $k_r$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('Decrement Distance Threshhold, m', 'Interpreter','latex', 'FontSize',14)
ylabel('Initial Call Duration, ms', 'Interpreter','latex', 'FontSize',14)
set(gca, 'TickLabelInterpreter', 'latex');
% set(gcf,'units','inch','position',[0,0,9,4]); 
figName = 'call_duration_threshold';
figpath = fullfile('/Users/ravi/Documents/projects/thesis/fig', figName);
% exportgraphics(gcf, [figpath '.png'], 'Resolution', 300, 'Append', false)
formatLatex(gca)

%% exporting to tikx format
% matlab2tikz('dist.tex', ...
%     'showLegend', false, 'width', '6.5cm', 'height', '4.2cm', ...
%     'extraAxisOptions', {'legend entries={none}', 'legend style={draw=none}', 'legend style={font=\footnotesize}'});
% stripLegendEntries('dist.tex');
%% model for call duration contraction due to doppler shift as a function of relative velocities

cd = 0.001:0.001:0.01;
vr = 1:50;
fs = 192e3;
c = 343;
for i = 1: length(cd)
    call = generateVirtualBatCall(20000, 90000, cd(i), fs, 20);
    for j = vr
        doppler_call = dopplerResample(call, -j, 343);
        call_contraction(i,j) = (length(call) - length(doppler_call))/fs;
    end
end
new_cd = 100 - 100*(cd' - call_contraction)./cd';
figure,
color_range = linspace(0, 1, length(cd));
for i = 1: length(cd)
    plot(vr, smooth(new_cd(i,:)), 'Color', [0, color_range(i),1-color_range(i)], 'LineWidth',0.75)
    hold on
end

grid on, grid minor

%%
axisbold
grid on
grid minor
ylim([0,15])
title('Signal Time Contraction due to Doppler Effect', 'Interpreter', 'latex');
xlabel('Relative Velocity, $v_r$, m/s', 'Interpreter','latex')
ylabel('Decrement in Original Duration, $\%$', 'Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'units','inch','position',[0,0,9,4]); 
figName = 'doppler_contraction';
figpath = fullfile('/Users/ravi/Documents/projects/thesis/fig', figName);
% exportgraphics(gcf, [figpath '.png'], 'Resolution', 300, 'Append', false)
%%

