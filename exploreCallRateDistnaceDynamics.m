%% Explore call rate versus distance across k_r values
% Minimal exploratory figure: C_r vs distance for selected k_r values.
% the script name may be changed to indicate the actual implementation
% Figures are displayed by default. Set savePlots = true after tuning.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'fcn'));
outDir = fullfile(thisDir, 'results', 'exploratory_figures');
savePlots = false;
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
distance_m = linspace(0.15, 6, 700);

curveColors = parula(numel(krValues));

%% Plot C_r versus distance
fig = figure('Color','w', 'Position', [160 160 500 500]);
hold on; box on; grid on;

CrReference_Hz = (c + vr) ./ (2 * (1 + 0) .* distance_m);
plot(distance_m, CrReference_Hz, ':', 'LineWidth', 1.5, ...
    'Color', [0.45 0.45 0.45], ...
    'DisplayName', '$k_r=0$');

for i = 1:numel(krValues)
    kr = krValues(i);
    Cr_Hz = (c + vr) ./ (2 * (1 + kr) .* distance_m);

    plot(distance_m, Cr_Hz, '-', 'LineWidth', 1.8, ...
        'Color', curveColors(i,:), ...
        'DisplayName', sprintf('$k_r=%d$', kr));
end

xlim([0 6]);
xlabel('Effective anchor distance, $d$ (m)', 'Interpreter','latex');
ylabel('Call rate, $C_r$ (Hz)', 'Interpreter','latex');
title('\textbf{Call-rate--distance relationship by $k_r$}', ...
    'Interpreter','latex');
axis square
grid minor
ylim([0 300]);
legend('Location','northeast', 'Interpreter','latex', 'FontSize', 9, 'Box','on');

eqText = { ...
    '$\displaystyle C_{r,\max}=\frac{c+v_r}{2(1+k_r)d_{\mathrm{stop}}}$', ...
    '$\displaystyle k_r=\frac{c+v_r}{2d_{\mathrm{stop}}C_{r,\max}}-1$'};
text(0.39, 0.52, eqText, ...
    'Units','normalized', ...
    'Interpreter','latex', ...
    'FontSize', 8, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');

formatLatex(fig, "half-square");

if savePlots
    exportPaperFigure(fig, fullfile(outDir, 'explore_Cr_vs_distance_by_kr'));
end
