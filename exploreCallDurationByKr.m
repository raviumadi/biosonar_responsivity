%% Explore call-duration evolution across k_r values
% Minimal exploratory figure: call duration dynamics under a no-overlap
% constraint, starting from a 5 ms call at maximum distance.
%
% Figures are displayed by default. Set savePlots = true after tuning.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'fcn'));
outDir = fullfile(thisDir, 'results', 'exploratory_figures');
savePlots = true;
runConfig = responsivityRunConfig();
if runConfig.OverrideOutputSwitches
    savePlots = runConfig.SaveFigures;
end
if runConfig.CloseFiguresBeforeRun
    close all;
end

if savePlots && ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Parameters
c = 343;                  % speed of sound, m/s
vr = 0;                   % stationary-anchor first-arrival reference
krValues = [3 5 7 9];
maxDistance_m = 5;
stopDistance_m = 0.15;
initialCallDuration_s = 0.003;
minCallDuration_s = 0.0001;

distance_m = linspace(stopDistance_m, maxDistance_m, 700);
curveColors = parula(numel(krValues));

%% Plot call-duration evolution
fig = figure('Color','w', 'Position', [160 160 500 500]);
hold on; box on; grid on;

for i = 1:numel(krValues)
    kr = krValues(i);

    Tdelay_s = 2 .* distance_m ./ (c + vr);
    Tcall_s = min(initialCallDuration_s, Tdelay_s ./ kr);
    Tcall_s = max(Tcall_s, minCallDuration_s);
    contractionDistance_m = (c + vr) * kr * initialCallDuration_s / 2;

    plot(distance_m, 1000 .* Tcall_s, '-', 'LineWidth', 1.8, ...
        'Color', curveColors(i,:), ...
        'DisplayName', sprintf('$k_r=%d$', kr));

    plot([contractionDistance_m contractionDistance_m], ...
        [0 1000 * initialCallDuration_s], '-.', ...
        'Color', curveColors(i,:), ...
        'LineWidth', 0.8, ...
        'HandleVisibility','off');
end

yline(1000 * initialCallDuration_s, '--', 'Color', [0.25 0.25 0.25], ...
    'LineWidth', 1.0, 'DisplayName', ...
    sprintf('$T_{c,0}=%.1f$ ms', 1000 * initialCallDuration_s));
yline(1000 * minCallDuration_s, '-.', 'Color', [0.45 0.45 0.45], ...
    'LineWidth', 1.0, 'DisplayName', ...
    sprintf('$T_{c,\\min}=%.1f$ ms', 1000 * minCallDuration_s));

xlim([0 maxDistance_m]);
ylim([0 6]);
xlabel('Effective anchor distance, $d$ (m)', 'Interpreter','latex', 'FontName', 'FixedWidth');
ylabel('Call duration, $T_c$ (ms)', 'Interpreter','latex');
title('\textbf{Call-duration contraction -- no-overlap constraint}', ...
    'Interpreter','latex');
axis square
grid minor
legend('Location','southeast', 'Interpreter','latex', 'FontSize', 8, 'Box','on');

eqText = {'$\displaystyle T^{\mathrm{delay}}=\frac{2d}{c+v_r}$', ...
    '$\displaystyle T_c=\max\left(T_{c,\min},\,\min\left[T_{c,0},\,\frac{T^{\mathrm{delay}}}{k_r}\right]\right)$', ...
    '$\displaystyle d_{\mathrm{contract}}=\frac{(c+v_r)k_rT_{c,0}}{2}$'};
text(0.08, 0.92, eqText, ...
    'Units','normalized', ...
    'Interpreter','latex', ...
    'FontSize', 9, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');

formatLatex(fig, "half-square");

if savePlots
    exportPaperFigure(fig, fullfile(outDir, 'explore_Tc_vs_distance_by_kr'));
end
