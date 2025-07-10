% MATLAB script to simulate distance covered before superposition
clear; clc; close all;

% Define chaser velocities (m/s)
chaser_speeds = 5;

% Define relative target velocities (m/s)
target_relative_speeds = linspace(0.1, 5, 500);

% Initial separation distance (m)
initial_distance = 10;

% Time step and max simulation time
dt = 0.01;  % seconds
max_time = 50; % seconds

% Initialize figure
figure; hold on;
colors = lines(length(target_relative_speeds));

for i = 1:length(chaser_speeds)
    chaser_speed = chaser_speeds(i);
    distances = zeros(1, length(target_relative_speeds));
    
    for j = 1:length(target_relative_speeds)
        target_speed = chaser_speed - target_relative_speeds(j);
        distance = initial_distance;
        time_elapsed = 0;
        
        % Simulate motion
        while distance > 0 && time_elapsed < max_time
            distance = distance - (chaser_speed - target_speed) * dt;
            time_elapsed = time_elapsed + dt;
        end
        
        distances(j) = chaser_speed * time_elapsed; % Distance covered by chaser
        timespan(j) = time_elapsed;
    end
    
    % Plot results
   time_plot_handle = plot(timespan, target_relative_speeds, '-.', "Color", 'k', 'LineWidth',2);
end

% Customize plot
ylabel('Relative Target Speed (m/s)', 'FontSize',14, Interpreter='latex');
xlabel('Time to Interception (s)','FontSize',14, Interpreter='latex');
title('Time to Interception for a Range of Relative Velocities', 'FontSize',14, Interpreter='latex');
subtitle('Bat Velocity = 5m/s. Initial Distance = 10m','FontSize', 12, Interpreter='latex')
% legend show;
ylim([.5*min(target_relative_speeds), max(target_relative_speeds)]);
xlim([min(timespan)-1 max(timespan)]);
% grid on;

% set(gcf,'units','normalized','position',[0.1300    0.1100    0.3    0.3]);
% grid minor
formatLatex(gca)
figpath = '~/Documents/projects/thesis/fig/chase_distance_model';
%%
exportgraphics(gcf, fullfile([figpath '.png']), 'Resolution', 300, 'BackgroundColor','none')
