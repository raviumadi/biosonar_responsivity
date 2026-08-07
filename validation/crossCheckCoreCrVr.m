%% Cross-check responsivity identities using simulateResponsivityCore
% Purpose: validate the internal consistency of the motion-aware
% responsivity identities using the new core simulator.
%
% This script checks the first-arrival special case:
%   phi = 0, no lag, no response floor, no call-rate ceiling
%
% Under those assumptions:
%   IPI = (1 + kr) * Tdelay
%   Tdelay = 2d / (c + vr)
% therefore:
%   C_r = (c + vr) / [2(1 + kr)d]

clear;
% clc;

thisDir = fileparts(mfilename('fullpath'));
simRoot = fileparts(thisDir);

addpath(simRoot);
addpath(fullfile(simRoot, 'fcn'));

paths = struct('figDir', fullfile(simRoot, 'results', 'validation_figures'));
savePlot = true;

% Choose validation mode:
%   "exact"        - deterministic first-arrival identity check; assertions on.
%   "jitteredBat"  - variable bat speed; useful for a more dynamic figure.
%   "movingTarget" - stochastic target motion; useful for stress-testing local updates.
validationScenario = "jitteredBat";
runConfig = responsivityRunConfig();
if runConfig.OverrideOutputSwitches
    savePlot = runConfig.SaveFigures;
    if runConfig.ValidationScenario ~= ""
        validationScenario = runConfig.ValidationScenario;
    end
end
if runConfig.CloseFiguresBeforeRun
    close all;
end
strictIdentityAssertions = validationScenario == "exact";

%% ------------ Run core simulation ------------
params = struct();
params.c = 343;
params.kr = 5;
params.initialDistance_m = 10;
params.initialBatSpeed_m_s = 8.0;
params.initialCallDuration_s = 0.005;
params.interceptDistance_m = 0.05;
params.maxCalls = 250;

opts = struct();
opts.geometryMode = "1D";
opts.numTargets = 1;
opts.targetMotion = false;
opts.batVelocityMode = "constant";
opts.anchorMode = "single";
opts.timingMode = "motionAware";
opts.phi = 0;
opts.lag.enabled = false;
opts.enforceMaxCallRate = false;
opts.enforceResponseFloor = false;
opts.callDurationMode = "constant";

switch validationScenario
    case "exact"
        % Leave deterministic settings unchanged.

    case "jitteredBat"
        opts.rngSeed = 10;
        opts.batVelocityMode = "jittered";
        opts.batVelocityJitterFrac = 0.20;

    case "movingTarget"
        opts.rngSeed = 11;
        opts.geometryMode = "1D";
        opts.targetMotion = true;
        opts.targetVelocityMode = "stochastic";
        opts.targetVelocityScale = 0.35;
        opts.targetVelocityJitterFrac = 0.20;

    otherwise
        error('Unknown validationScenario: %s', validationScenario);
end

res = simulateResponsivityCore(params, opts);
T = res.calls;

%% ------------ Extract arrays ------------
c = params.c;
kr = params.kr;

dt = T.IPI_s(:);
Cr = T.CallRate_Hz(:);
d_before = T.AnchorDistance_m(:);
Ta_step = T.Ta_s(:); %#ok<NASGU>
vr_before = T.RelativeClosingVelocityForTiming_m_s(:);
t_calls = T.CallOnsetTime_s(:);
t_step = T.NextCallOnsetTime_s(:);

%% ------------ Identities / theory ------------
% 1) Distance-from-call-rate identity in the motion-aware first-arrival form:
%    d = (c + v_r) / [2(1+k_r) C_r]
d_fromCr = (c + vr_before) ./ (2 * (1 + kr) .* Cr);

% 2) Relative velocity from the onset-to-onset step relation.
% The predicted displacement is the local-constant radial update used by the
% model; the observed displacement is the 3D position-derived distance change.
v_from_predicted_step = T.DistanceChangePredicted_m(:) ./ dt;
v_from_observed_step = T.DistanceChangeObserved_m(:) ./ dt;

% 3) Relative velocity from the motion-aware call-rate identity:
%    C_r = (c + v_r) / [2(1+k_r)d]
% => v_r = 2(1+k_r)d C_r - c
v_from_identity = 2 * (1 + kr) .* d_before .* Cr - c;

%% ------------ Plots ------------
fig = figure('Color','w', 'Position', [100 100 600 600]);

% (I) Relative velocity consistency over time
subplot(2,2,1);
plot(t_step, v_from_predicted_step, 'k-', 'LineWidth', 1.5); hold on;
plot(t_step, vr_before, 'r--', 'LineWidth', 1.5);
plot(t_step, v_from_observed_step, 'Color', [0.1 0.4 0.9], 'LineStyle', '-.', 'LineWidth', 1.5);
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Relative velocity $v_r$ (m/s)', 'Interpreter', 'latex');
lgd = legend('Predicted step $(\Delta d_{\rm pred}/\Delta t)$', ...
    'Stored timing $v_r$', 'Observed step $(\Delta d_{\rm obs}/\Delta t)$', ...
    'Location', 'best', 'Interpreter', 'latex');
lgd.Box = 'off';
title('$\textbf{I.}$ Relative velocity consistency','Interpreter', 'latex');
ylim([min(v_from_predicted_step)-0.5, max(v_from_predicted_step)+0.5]);
% formatLatex(gca);

% (II) Distance identity over time
subplot(2,2,2);
plot(t_calls, d_before, 'k-', 'LineWidth', 1.5); hold on;
plot(t_calls, d_fromCr, 'b--', 'LineWidth', 2);
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Distance to anchor $d$ (m)', 'Interpreter', 'latex');
lgd = legend('Simulated $d$ before step', '$d = (c+v_r)/[2(1+k_r)C_r]$', ...
    'Location', 'best', 'Interpreter', 'latex');
lgd.Box = 'off';
title('$\textbf{II.}$ Distance identity', 'Interpreter', 'latex');
% formatLatex(gca);

% (III) C_r versus d
subplot(2,2,3);
plot(d_before, Cr, 'k.', 'MarkerSize', 20); hold on;
plot(d_fromCr, Cr, 'r.', 'MarkerSize', 10);
xlabel('Distance $d$ (m)', 'Interpreter', 'latex');
ylabel('Call rate $C_r$ (Hz)', 'Interpreter', 'latex');
title('$\textbf{III.}$ $C_r(d) = (c+v_r) / [2(1+k_r)d]$', 'Interpreter', 'latex');
lgd = legend('Simulated pairs $(d,C_r)$', 'Identity using $C_r$', ...
    'Location', 'best', 'Interpreter', 'latex');
lgd.Box = 'off';
ylim([0 220]);
% formatLatex(gca);

% (IV) v_r versus C_r
subplot(2,2,4);
plot(Cr, v_from_predicted_step, 'k.', 'MarkerSize', 8); hold on;
plot(Cr, v_from_identity, 'ro', 'MarkerSize', 6);
xlabel('Call rate $C_r$ (Hz)', 'Interpreter', 'latex');
ylabel('Relative velocity $v_r$ (m/s)', 'Interpreter', 'latex');
title('$\textbf{IV.}$ $v_r$ vs $C_r$', 'Interpreter', 'latex');
lgd = legend('From core step', 'From motion-aware identity', ...
    'Location', 'best', 'Interpreter', 'latex');
lgd.Box = 'off';
ylim([min(v_from_predicted_step)-0.5, max(v_from_predicted_step)+0.5]);
xlim([0 220]);
% formatLatex(gca);

axesHandles = findall(fig, 'Type', 'axes');
for ax = reshape(axesHandles, 1, [])
    grid(ax, 'on');
    grid(ax, 'minor');
    box(ax, 'on');
end

mainTitle = sgtitle('Responsivity core: validation of derived parameters', ...
    'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'latex');

paperStyle = formatLatex(fig, "compact-square");
% mainTitle.FontSize = paperStyle.FigureTitleFontSize;
% mainTitle.FontWeight = 'bold';

%% ------------ Sanity metrics ------------
meanVrErr = mean(abs(vr_before - v_from_predicted_step), 'omitnan');
medianDistErr = median(abs(d_before - d_fromCr), 'omitnan');
medianVrIdentityErr = median(abs(v_from_predicted_step - v_from_identity), 'omitnan');
medianObservedStepErr = median(abs(v_from_observed_step - v_from_predicted_step), 'omitnan');

fprintf('Mean |v_r(profile) - v_r(predicted step)|: %.6e m/s\n', meanVrErr);
fprintf('Median |d_before - d_fromCr|: %.6e m\n', medianDistErr);
fprintf('Median |v_from_step - v_from_identity|: %.6e m/s\n', medianVrIdentityErr);
fprintf('Median |v_observed_step - v_predicted_step|: %.6e m/s\n', medianObservedStepErr);

if strictIdentityAssertions
    assert(meanVrErr < 1e-10, 'Stored timing velocity and predicted step velocity diverged.');
    assert(medianDistErr < 1e-10, 'Distance identity failed for the first-arrival core case.');
    assert(medianVrIdentityErr < 1e-10, 'Velocity identity failed for the first-arrival core case.');
    disp('Core responsivity identity cross-check passed.');
else
    disp('Core responsivity exploratory cross-check completed.');
    disp('Note: hard identity assertions are disabled outside the exact deterministic scenario.');
end

%% ------------ Save figure ------------
if savePlot
    if ~exist(paths.figDir, 'dir')
        mkdir(paths.figDir);
    end
    figPath = char(fullfile(paths.figDir, "core_equation_validation_" + validationScenario));
    exportPaperFigure(fig, figPath);
end
