% delta plots - relationships between sonar parameters and motion
% 1. delta s as a function of subject velocity and target distance
c = 343;
target_distance = linspace(0.02,2, 200);
subject_velocity = linspace(1, 10, 10);

% Initialize arrays to store results
deltaS_results = zeros(length(subject_velocity), length(target_distance));
signal_flight_time_results = zeros(length(subject_velocity), ...
    length(target_distance));
new_target_distance_results = zeros(length(subject_velocity), ...
    length(target_distance));

% Loop over subject_velocity
for i = 1:length(subject_velocity)
    % Loop over target_distance
    for j = 1:length(target_distance)
        % Calculate deltaS, signal_flight_time, and new_target_distance
        [deltaS, signal_flight_time, new_target_distance] = ...
            emissionReceptionSpaceTimeChange(subject_velocity(i), ...
            target_distance(j));
        
        % Store results in arrays
        deltaS_results(i, j) = deltaS;
        signal_flight_time_results(i, j) = signal_flight_time;
        new_target_distance_results(i, j) = new_target_distance;
    end
end

% calculate call rate for different Ta:Tb ratios- CR = ((Ta+Tb)/2) s
CR = 1./signal_flight_time_results;
CR_quarter = 1./(signal_flight_time_results +signal_flight_time_results./4 );
CR_half = 1./(signal_flight_time_results + signal_flight_time_results./2);
CR_equal = 1./(2*signal_flight_time_results);
CR_double = 1./(signal_flight_time_results + 2 * signal_flight_time_results);
CR_triple = 1./(signal_flight_time_results + 3 * signal_flight_time_results);
CR_quad = 1./(signal_flight_time_results + 4 * signal_flight_time_results);
CR_penta = 1./(signal_flight_time_results + 5 * signal_flight_time_results);
CR_hexa = 1./(signal_flight_time_results + 6 * signal_flight_time_results);
CR_hepta = 1./(signal_flight_time_results + 7 * signal_flight_time_results);
CR_octa = 1./(signal_flight_time_results + 8 * signal_flight_time_results);
CR_nona = 1./(signal_flight_time_results + 9 * signal_flight_time_results);
CR_deca = 1./(signal_flight_time_results + 10 * signal_flight_time_results);
% plot
plot(target_distance, CR, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', 'k')
hold on
plot(target_distance, CR_quarter, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', '#f01c05')
plot(target_distance, CR_half, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', '#f05b05')

plot(target_distance, CR_equal, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', '#f0a505')
plot(target_distance, CR_double, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', '#9af005')
plot(target_distance, CR_triple, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', '#05f038')
plot(target_distance, CR_quad, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', '#05f09d')
plot(target_distance, CR_penta, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', '#05e4f0')
plot(target_distance, CR_hexa, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', '#055ff0')
plot(target_distance, CR_hepta, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', '#1405f0')
plot(target_distance, CR_octa, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', '#5b05f0')
plot(target_distance, CR_nona, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', '#0543f0')
plot(target_distance, CR_deca, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', '#0d039c')
grid on, 
% set(gca, 'TickLabelInterpreter', 'latex');
ylim([0 500])
% axisbold

%
% legend('0', '0.25', '0.5', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10')
title('Call Rate as a Function of $T_a$ \& $T_b$', 'Interpreter','latex', 'FontSize',14)
xlabel('Target Distance, m', 'Interpreter','latex', 'FontSize',14)
ylabel('Calls/s', 'Interpreter','latex', 'FontSize',14)
formatLatex(gca)
grid minor
% set(gcf,'units','inch','position',[0,0,9,4]); 
figName = 'call_rate_limit';
figpath = fullfile('/Users/ravi/Documents/projects/thesis/fig', figName);
% exportgraphics(gcf, [figpath '.png'], 'Resolution', 300, 'Append', false)

formatLatex(gca)
% 
% % set(fig,'visible','off')
% matlab2tikz('cr.tex', ...
%     'showLegend', false, 'width', '6.5cm', 'height', '4.2cm', ...
%     'extraAxisOptions', {'legend entries={none}', 'legend style={draw=none}', 'legend style={font=\footnotesize}'});
% stripLegendEntries('cr.tex');
%%
% Plotting deltaS against target_distance for each subject_velocity
figure;
hold on;
for i = 1:length(subject_velocity)
    plot(target_distance, deltaS_results(i, :)./343, 'DisplayName', ...
        sprintf('Subject Velocity = %.2f m/s', subject_velocity(i)));
end
hold off;
xlabel('Target Distance (m)');
ylabel('DeltaS (m)');
title(['DeltaS as a Function of Target Distance for' ...
    ' Different Subject Velocities']);
% legend('Location', 'best');
% grid minor

%%
figure;
hold on;
for i = 1:size(signal_flight_time_results,1)
    plot(target_distance, signal_flight_time_results(i, :));
end
hold off;
xlabel('Target Distance (m)');
ylabel('Signal Return Time, s');
title(['Signal Return Time as a Function of Target Distance for' ...
    ' Different Subject Velocities']);
grid minor

