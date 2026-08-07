%% Explore buzz readiness from the responsivity curve
% This script uses the responsivity-curve definition already set out in the
% manuscript:
%   R_n = 1 / |Delta t_(n+1) - Delta t_n|
% and defines buzz readiness as the call index at which R_n comes closest
% to the call-rate ceiling C_r,max.
%
% The sweep is designed to ask how the inferred buzz-readiness distance
% changes with:
%   1. responsivity coefficient k_r,
%   2. starting velocity v_0,
%   3. starting call duration T_c,0.
%
% The simulations use a one-target, three-dimensional, motion-aware
% setting with a stationary target and a non-jittered bat so the resulting
% trajectories are easier to interpret.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir, 'fcn'));

outDir = fullfile(thisDir, 'results', 'exploratory_figures');
savePlots = false;
saveStats = false;
runConfig = responsivityRunConfig();
if runConfig.OverrideOutputSwitches
    savePlots = runConfig.SaveFigures;
    saveStats = runConfig.SaveTables;
end
if runConfig.CloseFiguresBeforeRun
    close all;
end
showExampleCurves = true;
ceilingHitFrac = 0.99;
terminalProxyMode = "TbMin"; % "TbMin" or "IpiMin"

if (savePlots || saveStats) && ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Sweep settings
krValues = [3 5 7 9];
startVelocityValues_m_s = 1:5;
initialCallDurationValues_ms = [3 5 7 10];
nRepeatsPerCondition = 20;

dStart_m = 5;
dStop_m = 0.15;
phiFixed = 0.1;

%% Shared simulation settings
paramsBase = struct();
paramsBase.c = 343;
paramsBase.kr = 5;
paramsBase.initialDistance_m = dStart_m;
paramsBase.initialBatSpeed_m_s = 3;
paramsBase.initialCallDuration_s = 0.005;
paramsBase.minCallDuration_s = 0.0005;
paramsBase.maxCalls = 300;
paramsBase.maxElapsedTime_s = 10;
paramsBase.maxAnchorDistance_m = 10;
paramsBase.interceptDistance_m = dStop_m;
paramsBase.numSequences = 1;

optsBase = struct();
optsBase.geometryMode = "3D";
optsBase.numTargets = 1;
optsBase.targetMotion = false;
optsBase.targetVelocityMode = "constant";
optsBase.targetVelocityScale = 0;
optsBase.targetVelocityJitterFrac = 0;
optsBase.batVelocityMode = "constant";
optsBase.batVelocityJitterFrac = 0;
optsBase.anchorMode = "single";
optsBase.timingMode = "motionAware";
optsBase.callDurationMode = "previousTa";
optsBase.enforceMaxCallRate = true;
optsBase.callDurationJitter.enabled = false;
optsBase.callDurationJitter.mode = "additive";
optsBase.callDurationJitter.rho = 0.25;
optsBase.phi = phiFixed;

rng(42);

%% Run sweep
sequenceSummary = table();
curveTable = table();
seqCounter = 0;

for iKr = 1:numel(krValues)
    kr = krValues(iKr);

    for iVel = 1:numel(startVelocityValues_m_s)
        v0 = startVelocityValues_m_s(iVel);

        for iTc0 = 1:numel(initialCallDurationValues_ms)
            tc0_ms = initialCallDurationValues_ms(iTc0);

            for r = 1:nRepeatsPerCondition
                seqCounter = seqCounter + 1;

                params = paramsBase;
                params.kr = kr;
                params.initialBatSpeed_m_s = v0;
                params.initialCallDuration_s = tc0_ms / 1000;
                params.maxCallRate_Hz = deriveCrMaxFromStop(params.c, kr, dStop_m);

                opts = optsBase;
                opts.rngSeed = [];

                res = simulateResponsivityCore(params, opts);
                T = res.calls;
                if isempty(T) || height(T) < 3
                    continue
                end

                [seqRow, curveRows] = computeResponsivityMetrics(T, seqCounter, kr, v0, tc0_ms, params.maxCallRate_Hz, dStop_m, ceilingHitFrac);
                sequenceSummary = appendCompatible(sequenceSummary, seqRow);
                curveTable = appendCompatible(curveTable, curveRows);
            end
        end
    end
end

fprintf('Total sequences analysed: %d\n', height(sequenceSummary));
fprintf('Responsivity-curve points retained: %d\n', height(curveTable));
fprintf('k_r values: %s\n', strjoin(string(krValues), ', '));
fprintf('Starting velocities: %s m/s\n', strjoin(string(startVelocityValues_m_s), ', '));
fprintf('Initial call durations: %s ms\n', strjoin(string(initialCallDurationValues_ms), ', '));
fprintf('Call-rate ceilings derived from d_stop = %.2f m\n', dStop_m);

%% Summaries
summaryByKr = summariseReadiness(sequenceSummary, "kr", krValues);
summaryByVelocity = summariseReadiness(sequenceSummary, "InitialVelocity_m_s", startVelocityValues_m_s);
summaryByDuration = summariseReadiness(sequenceSummary, "InitialCallDuration_ms", initialCallDurationValues_ms);

fprintf('\n=== Buzz-readiness summary by k_r ===\n');
disp(summaryByKr);

fprintf('\n=== Buzz-readiness summary by starting velocity ===\n');
disp(summaryByVelocity);

fprintf('\n=== Buzz-readiness summary by initial call duration ===\n');
disp(summaryByDuration);

ttcCompareByKr = compareTimeToContactProxy(sequenceSummary, "kr", krValues);
ttcCompareByVelocity = compareTimeToContactProxy(sequenceSummary, "InitialVelocity_m_s", startVelocityValues_m_s);

proxyConfig = getTerminalProxyConfig(terminalProxyMode);
terminalProxyByKr = compareTerminalProxy(sequenceSummary, "kr", krValues, proxyConfig);
terminalProxyByVelocity = compareTerminalProxy(sequenceSummary, "InitialVelocity_m_s", startVelocityValues_m_s, proxyConfig);

fprintf('\n=== Paired comparison: \\delta t* vs %s by k_r ===\n', proxyConfig.PrintLabel);
disp(terminalProxyByKr);

fprintf('\n=== Paired comparison: \\delta t* vs %s by starting velocity ===\n', proxyConfig.PrintLabel);
disp(terminalProxyByVelocity);

fprintf('\n=== Paired comparison: time to target from buzz readiness vs from d_stop by k_r ===\n');
disp(ttcCompareByKr);

fprintf('\n=== Paired comparison: time to target from buzz readiness vs from d_stop by starting velocity ===\n');
disp(ttcCompareByVelocity);

if saveStats
    writetable(sequenceSummary, fullfile(outDir, 'responsivity_curve_sequence_metrics.csv'));
    writetable(curveTable, fullfile(outDir, 'responsivity_curve_all_points.csv'));
    writetable(summaryByKr, fullfile(outDir, 'responsivity_curve_summary_by_kr.csv'));
    writetable(summaryByVelocity, fullfile(outDir, 'responsivity_curve_summary_by_velocity.csv'));
    writetable(summaryByDuration, fullfile(outDir, 'responsivity_curve_summary_by_duration.csv'));
    writetable(terminalProxyByKr, fullfile(outDir, sprintf('responsivity_curve_%s_proxy_by_kr.csv', proxyConfig.FileStem)));
    writetable(terminalProxyByVelocity, fullfile(outDir, sprintf('responsivity_curve_%s_proxy_by_velocity.csv', proxyConfig.FileStem)));
    writetable(ttcCompareByKr, fullfile(outDir, 'responsivity_curve_ttc_proxy_by_kr.csv'));
    writetable(ttcCompareByVelocity, fullfile(outDir, 'responsivity_curve_ttc_proxy_by_velocity.csv'));
end

%% Main figure: buzz-readiness distance summaries
colorsKr = [
    0.28 0.18 0.72;
    0.18 0.58 0.90;
    0.41 0.76 0.31;
    0.94 0.86 0.10
];

fig = figure('Color', 'w', 'Position', [80 90 1000 1000]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, '\textbf{Buzz readiness from the responsivity curve in single-target pursuit simulations}', ...
    'Interpreter', 'latex', 'FontSize', 15);

nexttile;
plotReadinessByKr(sequenceSummary, krValues, colorsKr);
title({'\textbf{A. Responsivity coefficient}', ...
    sprintf('$d_0=%.1f$ m, $\\phi=%.1f$, $v_0=1$--$5$ m s$^{-1}$, $T_{c,0}=3$--$10$ ms', dStart_m, phiFixed)}, ...
    'Interpreter', 'latex');
axis square

nexttile;
plotReadinessTrend(sequenceSummary, "InitialVelocity_m_s", startVelocityValues_m_s, ...
    "Median buzz-readiness distance (m)", colorsKr, "kr", krValues, '$v_0$ (m s$^{-1}$)');
title({'\textbf{B. Starting velocity}', ...
    sprintf('$d_0=%.1f$ m, $\\phi=%.1f$, $T_{c,0}=3$--$10$ ms', dStart_m, phiFixed)}, ...
    'Interpreter', 'latex');
axis square

nexttile;
plotReadinessPercentViolin(sequenceSummary, krValues, startVelocityValues_m_s, colorsKr);
title({'\textbf{C. Buzz-readiness position within the sequence}', ...
    sprintf('$d_0=%.1f$ m, $\\phi=%.1f$, readiness call as \\%% of full seq.', dStart_m, phiFixed)}, ...
    'Interpreter', 'latex');
axis square

nexttile;
plotDeltaTStarViolin(sequenceSummary, krValues, startVelocityValues_m_s, colorsKr, ceilingHitFrac, proxyConfig);
title({'\textbf{D. Buzz-readiness timing interval, $\delta t^*$}', ...
    sprintf('$d_0=%.1f$ m, $\\phi=%.1f$, with terminal %s overlay (grey)', dStart_m, phiFixed, proxyConfig.TitleLabel)}, ...
    'Interpreter', 'latex');
axis square

formatLatex(fig, "full-square");

if savePlots
    exportPaperFigure(fig, fullfile(outDir, 'explore_responsivity_curve_buzz_readiness'));
end

%% Optional diagnostic figure: example responsivity curves
if showExampleCurves
    figCurve = figure('Color', 'w', 'Position', [100 120 1000 1000]);
    tl2 = tiledlayout(figCurve, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl2, '\textbf{Responsivity curves and time to contact at buzz readiness}', ...
        'Interpreter', 'latex', 'FontSize', 15);


    makeExampleCurvePanel(curveTable, sequenceSummary, "kr", krValues, colorsKr, ...
        fixedValueMask(sequenceSummary, "InitialVelocity_m_s", 3) & ...
        fixedValueMask(sequenceSummary, "InitialCallDuration_ms", 5), ...
        '\textbf{A. Varying $k_r$}', '$v_0=3$ m s$^{-1}$, $T_{c,0}=5$ ms');
    axis square
    ylim([0 10]);

    makeExampleCurvePanel(curveTable, sequenceSummary, "InitialVelocity_m_s", startVelocityValues_m_s, turbo(numel(startVelocityValues_m_s)), ...
        fixedValueMask(sequenceSummary, "kr", 5) & ...
        fixedValueMask(sequenceSummary, "InitialCallDuration_ms", 5), ...
        '\textbf{B. Varying $v_0$}', '$k_r=5$, $T_{c,0}=5$ ms');
    axis square
    ylim([0 10]);

    makeExampleCurvePanel(curveTable, sequenceSummary, "InitialCallDuration_ms", initialCallDurationValues_ms, parula(numel(initialCallDurationValues_ms)), ...
        fixedValueMask(sequenceSummary, "kr", 5) & ...
        fixedValueMask(sequenceSummary, "InitialVelocity_m_s", 3), ...
        '\textbf{C. Varying $T_{c,0}$}', '$k_r=5$, $v_0=3$ m s$^{-1}$');
    axis square

    nexttile;
    plotReadinessTimeToContactViolin(sequenceSummary, krValues, startVelocityValues_m_s, colorsKr);
    title({'\textbf{D. Time to contact}', ...
        sprintf('$d_0=%.1f$ m, $\\phi=%.1f$, buzz readiness to target with $d_{\\mathrm{stop}}$ reference', dStart_m, phiFixed)}, ...
        'Interpreter', 'latex');
    axis square

    formatLatex(figCurve, "full-square");

    if savePlots
        exportPaperFigure(figCurve, fullfile(outDir, 'explore_responsivity_curve_examples'));
    end
end

%% Local functions
function crMax = deriveCrMaxFromStop(c, kr, dStop)
    crMax = c ./ (2 * (1 + kr) * dStop);
end

function [seqRow, curveRows] = computeResponsivityMetrics(T, seqID, kr, v0, tc0_ms, crMax, dStop_m, ceilingHitFrac)
    ipi = T.IPI_s(:);
    d = T.AnchorDistance_m(:);
    callRate = T.CallRate_Hz(:);

    deltaIpi = abs(diff(ipi));
    valid = isfinite(deltaIpi) & deltaIpi > 0 & ...
        isfinite(d(1:end-1)) & isfinite(d(2:end)) & ...
        isfinite(callRate(1:end-1));

    resp = nan(size(deltaIpi));
    resp(valid) = 1 ./ deltaIpi(valid);
    respDistance = nan(size(deltaIpi));
    respDistance(valid) = 0.5 * (d(1:end-1) + d(2:end));
    deltaTstar_ms = 1000 * deltaIpi;

    score = abs(resp - crMax);
    score(~valid) = nan;
    [~, idxStar] = min(score);

    if isempty(idxStar) || ~isfinite(score(idxStar))
        seqRow = table();
        curveRows = table();
        return
    end

    readinessSequencePercent = 100 * T.CallNumber(idxStar) / max(height(T), 1);
    readinessCallIdx = max(1, min(height(T), T.CallNumber(idxStar)));
    tbTerminalWindow = T.TbEffective_s(readinessCallIdx:end);
    terminalTbMin_ms = 1000 * min(tbTerminalWindow(isfinite(tbTerminalWindow)), [], 'omitnan');
    ipiTerminalWindow = T.IPI_s(readinessCallIdx:end);
    terminalIpiMin_ms = 1000 * min(ipiTerminalWindow(isfinite(ipiTerminalWindow)), [], 'omitnan');
    maxObservedCallRate = max(callRate, [], 'omitnan');
    achievedCeiling = maxObservedCallRate >= ceilingHitFrac * crMax;
    if ismember('RelativeClosingVelocity_m_s', T.Properties.VariableNames)
        vrReadiness = T.RelativeClosingVelocity_m_s(readinessCallIdx);
    else
        vrReadiness = v0;
    end
    if ~isfinite(vrReadiness) || vrReadiness <= 0
        vrReadiness = v0;
    end
    timeToTargetAtReadiness_ms = 1000 * max(respDistance(idxStar), 0) / max(vrReadiness, eps);
    timeFromDstopToTarget_ms = 1000 * max(dStop_m, 0) / max(vrReadiness, eps);

    seqRow = table(seqID, kr, v0, tc0_ms, crMax, ...
        respDistance(idxStar), resp(idxStar), resp(idxStar) / crMax, ...
        deltaTstar_ms(idxStar), T.CallNumber(idxStar), readinessSequencePercent, ...
        terminalTbMin_ms, terminalIpiMin_ms, ...
        timeToTargetAtReadiness_ms, timeFromDstopToTarget_ms, achievedCeiling, maxObservedCallRate, height(T), ...
        'VariableNames', {'SeqID','kr','InitialVelocity_m_s','InitialCallDuration_ms', ...
        'CrMax_Hz','BuzzReadinessDistance_m','PrimeResponsivity_Hz', ...
        'PrimeResponsivityNorm','DeltaTStar_ms','BuzzReadinessCallNumber', ...
        'BuzzReadinessSequencePercent','TerminalTbMin_ms','TerminalIpiMin_ms', ...
        'TimeToTargetAtReadiness_ms','TimeFromDstopToTarget_ms','AchievedCeiling', ...
        'MaxObservedCallRate_Hz','NumCalls'});

    nRows = numel(resp);
    curveRows = table( ...
        repmat(seqID, nRows, 1), ...
        repmat(kr, nRows, 1), ...
        repmat(v0, nRows, 1), ...
        repmat(tc0_ms, nRows, 1), ...
        repmat(crMax, nRows, 1), ...
        (1:nRows)', ...
        respDistance, ...
        resp, ...
        resp / crMax, ...
        deltaTstar_ms, ...
        score, ...
        repmat(idxStar, nRows, 1), ...
        'VariableNames', {'SeqID','kr','InitialVelocity_m_s','InitialCallDuration_ms', ...
        'CrMax_Hz','ResponsivityIndex','ResponsivityDistance_m', ...
        'Responsivity_Hz','ResponsivityNorm','DeltaIpi_ms', ...
        'DistanceToCeiling_Hz','BuzzReadinessIndex'});
end

function summaryTbl = summariseReadiness(T, groupVar, groupValues)
    rowsOut = cell(numel(groupValues), 1);

    for i = 1:numel(groupValues)
        g = groupValues(i);
        rows = abs(T.(groupVar) - g) < 1e-12 & isfinite(T.BuzzReadinessDistance_m);
        S = T(rows, :);
        if isempty(S)
            rowsOut{i} = table(g, 0, nan, nan, nan, nan, nan, nan, nan, nan, nan, ...
                'VariableNames', {char(groupVar),'NumSequences','MedianBuzzReadinessDistance_m', ...
                'Q1BuzzReadinessDistance_m','Q3BuzzReadinessDistance_m', ...
                'MeanBuzzReadinessDistance_m','MedianDeltaTStar_ms','MeanDeltaTStar_ms', ...
                'MedianPrimeResponsivityNorm','MedianMaxObservedCallRate_Hz','MedianNumCalls'});
            continue
        end

        rowsOut{i} = table(g, height(S), ...
            median(S.BuzzReadinessDistance_m, 'omitnan'), ...
            prctile(S.BuzzReadinessDistance_m, 25), ...
            prctile(S.BuzzReadinessDistance_m, 75), ...
            mean(S.BuzzReadinessDistance_m, 'omitnan'), ...
            median(S.DeltaTStar_ms, 'omitnan'), ...
            mean(S.DeltaTStar_ms, 'omitnan'), ...
            median(S.PrimeResponsivityNorm, 'omitnan'), ...
            median(S.MaxObservedCallRate_Hz, 'omitnan'), ...
            median(S.NumCalls, 'omitnan'), ...
            'VariableNames', {char(groupVar),'NumSequences','MedianBuzzReadinessDistance_m', ...
            'Q1BuzzReadinessDistance_m','Q3BuzzReadinessDistance_m', ...
            'MeanBuzzReadinessDistance_m','MedianDeltaTStar_ms','MeanDeltaTStar_ms', ...
            'MedianPrimeResponsivityNorm','MedianMaxObservedCallRate_Hz','MedianNumCalls'});
    end

    summaryTbl = vertcat(rowsOut{:});
end

function proxyConfig = getTerminalProxyConfig(mode)
    switch string(mode)
        case "TbMin"
            proxyConfig = struct( ...
                'VariableName', "TerminalTbMin_ms", ...
                'PrintLabel', 'terminal T_b,min', ...
                'TitleLabel', '$T_{b,\min}$', ...
                'AxisLabel', 'terminal $T_{b,\min}$', ...
                'LegendLabel', 'minimum terminal $T_b$ reached after buzz readiness', ...
                'FileStem', 'tb');
        case "IpiMin"
            proxyConfig = struct( ...
                'VariableName', "TerminalIpiMin_ms", ...
                'PrintLabel', 'terminal IPI_min', ...
                'TitleLabel', '$\mathrm{IPI}_{\min}$', ...
                'AxisLabel', 'terminal $\mathrm{IPI}_{\min}$', ...
                'LegendLabel', 'minimum terminal IPI reached after buzz readiness', ...
                'FileStem', 'ipimin');
        otherwise
            error('Unknown terminalProxyMode: %s', string(mode));
    end
end

function summaryTbl = compareTerminalProxy(T, groupVar, groupValues, proxyConfig)
    rowsOut = cell(numel(groupValues), 1);

    for i = 1:numel(groupValues)
        g = groupValues(i);
        rows = abs(T.(groupVar) - g) < 1e-12 & ...
            isfinite(T.DeltaTStar_ms) & isfinite(T.(proxyConfig.VariableName));
        S = T(rows, :);

        if isempty(S)
            rowsOut{i} = table(g, 0, nan, nan, nan, nan, nan, nan, nan, nan, ...
                'VariableNames', {char(groupVar), 'NumSequences', ...
                'MedianDeltaTStar_ms', 'MedianTerminalProxy_ms', ...
                'MedianDifference_ms', 'MeanDifference_ms', ...
                'MADifference_ms', 'SignrankPValue', 'CohensDz', ...
                'PercentWithin10Percent'});
            continue
        end

        proxyVals = S.(proxyConfig.VariableName);
        delta = S.DeltaTStar_ms - proxyVals;
        denom = std(delta, 'omitnan');
        if ~isfinite(denom) || denom == 0
            cohensDz = nan;
        else
            cohensDz = mean(delta, 'omitnan') / denom;
        end

        relDiff = abs(delta) ./ max(proxyVals, eps);
        within10 = 100 * mean(relDiff <= 0.10, 'omitnan');

        try
            pSignrank = signrank(S.DeltaTStar_ms, proxyVals);
        catch
            pSignrank = nan;
        end

        rowsOut{i} = table(g, height(S), ...
            median(S.DeltaTStar_ms, 'omitnan'), ...
            median(proxyVals, 'omitnan'), ...
            median(delta, 'omitnan'), ...
            mean(delta, 'omitnan'), ...
            mean(abs(delta), 'omitnan'), ...
            pSignrank, cohensDz, within10, ...
            'VariableNames', {char(groupVar), 'NumSequences', ...
            'MedianDeltaTStar_ms', 'MedianTerminalProxy_ms', ...
            'MedianDifference_ms', 'MeanDifference_ms', ...
            'MADifference_ms', 'SignrankPValue', 'CohensDz', ...
            'PercentWithin10Percent'});
    end

    summaryTbl = vertcat(rowsOut{:});
end

function summaryTbl = compareTimeToContactProxy(T, groupVar, groupValues)
    rowsOut = cell(numel(groupValues), 1);

    for i = 1:numel(groupValues)
        g = groupValues(i);
        rows = abs(T.(groupVar) - g) < 1e-12 & ...
            isfinite(T.TimeToTargetAtReadiness_ms) & isfinite(T.TimeFromDstopToTarget_ms);
        S = T(rows, :);

        if isempty(S)
            rowsOut{i} = table(g, 0, nan, nan, nan, nan, nan, nan, nan, nan, ...
                'VariableNames', {char(groupVar), 'NumSequences', ...
                'MedianReadinessTTC_ms', 'MedianDstopTTC_ms', ...
                'MedianDifference_ms', 'MeanDifference_ms', ...
                'MADifference_ms', 'SignrankPValue', 'CohensDz', ...
                'PercentWithin10Percent'});
            continue
        end

        delta = S.TimeToTargetAtReadiness_ms - S.TimeFromDstopToTarget_ms;
        denom = std(delta, 'omitnan');
        if ~isfinite(denom) || denom == 0
            cohensDz = nan;
        else
            cohensDz = mean(delta, 'omitnan') / denom;
        end

        relDiff = abs(delta) ./ max(S.TimeFromDstopToTarget_ms, eps);
        within10 = 100 * mean(relDiff <= 0.10, 'omitnan');

        try
            pSignrank = signrank(S.TimeToTargetAtReadiness_ms, S.TimeFromDstopToTarget_ms);
        catch
            pSignrank = nan;
        end

        rowsOut{i} = table(g, height(S), ...
            median(S.TimeToTargetAtReadiness_ms, 'omitnan'), ...
            median(S.TimeFromDstopToTarget_ms, 'omitnan'), ...
            median(delta, 'omitnan'), ...
            mean(delta, 'omitnan'), ...
            mean(abs(delta), 'omitnan'), ...
            pSignrank, cohensDz, within10, ...
            'VariableNames', {char(groupVar), 'NumSequences', ...
            'MedianReadinessTTC_ms', 'MedianDstopTTC_ms', ...
            'MedianDifference_ms', 'MeanDifference_ms', ...
            'MADifference_ms', 'SignrankPValue', 'CohensDz', ...
            'PercentWithin10Percent'});
    end

    summaryTbl = vertcat(rowsOut{:});
end

function plotReadinessByKr(T, krValues, colorsKr)
    hold on; box on; grid on; grid minor;

    for i = 1:numel(krValues)
        kr = krValues(i);
        rows = abs(T.kr - kr) < 1e-12 & isfinite(T.BuzzReadinessDistance_m);
        y = filterRobustOutliers(T.BuzzReadinessDistance_m(rows));
        drawSimpleViolin(kr, y, colorsKr(i,:), 0.20);

        if ~isempty(y)
            q1 = prctile(y, 25);
            q3 = prctile(y, 75);
            med = median(y);
            plot([kr kr], [q1 q3], '-', 'Color', colorsKr(i,:), 'LineWidth', 2.2);
            plot(kr, med, 'o', 'MarkerFaceColor', colorsKr(i,:), ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 6);
        end
    end

    xlim([min(krValues)-0.7 max(krValues)+0.7]);
    xlabel('Responsivity coefficient, $k_r$', 'Interpreter', 'latex');
    ylabel('Buzz-readiness distance (m)', 'Interpreter', 'latex');
end

function plotReadinessTrend(T, groupVar, groupValues, yLabelText, colorsKr, lineVar, lineValues, xLabelText)
    hold on; box on; grid on; grid minor;

    for i = 1:numel(lineValues)
        lv = lineValues(i);
        xVals = nan(numel(groupValues), 1);
        medVals = nan(numel(groupValues), 1);
        q1Vals = nan(numel(groupValues), 1);
        q3Vals = nan(numel(groupValues), 1);

        for j = 1:numel(groupValues)
            gv = groupValues(j);
            rows = abs(T.(lineVar) - lv) < 1e-12 & ...
                abs(T.(groupVar) - gv) < 1e-12 & ...
                isfinite(T.BuzzReadinessDistance_m);
            y = T.BuzzReadinessDistance_m(rows);
            if isempty(y)
                continue
            end

            xVals(j) = gv;
            medVals(j) = median(y);
            q1Vals(j) = prctile(y, 25);
            q3Vals(j) = prctile(y, 75);
        end

        good = isfinite(xVals) & isfinite(medVals);
        if nnz(good) < 2
            continue
        end

        fill([xVals(good); flipud(xVals(good))], [q1Vals(good); flipud(q3Vals(good))], ...
            colorsKr(i,:), 'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(xVals(good), medVals(good), '-', 'Color', colorsKr(i,:), ...
            'LineWidth', 2.0, 'DisplayName', sprintf('$k_r=%g$', lv));
        scatter(xVals(good), medVals(good), 20, colorsKr(i,:), 'filled', 'HandleVisibility', 'off');
    end

    xlabel(xLabelText, 'Interpreter', 'latex');
    ylabel(yLabelText, 'Interpreter', 'latex');
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 8, 'Box', 'on');
end

function plotReadinessPercentViolin(T, krValues, velocityValues, colorsKr)
    hold on; box on; grid on; grid minor;

    xKr = 1:numel(krValues);
    xVel = (numel(krValues) + 2):(numel(krValues) + 1 + numel(velocityValues));
    violinWidth = 0.34;

    for i = 1:numel(krValues)
        rows = abs(T.kr - krValues(i)) < 1e-12 & isfinite(T.BuzzReadinessSequencePercent);
        y = T.BuzzReadinessSequencePercent(rows);
        drawSimpleViolin(xKr(i), y, colorsKr(i,:), violinWidth);
        if ~isempty(y)
            plot(xKr(i), median(y, 'omitnan'), 'o', 'MarkerFaceColor', colorsKr(i,:), ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');
        end
    end

    velColors = turbo(numel(velocityValues));
    for i = 1:numel(velocityValues)
        rows = abs(T.InitialVelocity_m_s - velocityValues(i)) < 1e-12 & isfinite(T.BuzzReadinessSequencePercent);
        y = T.BuzzReadinessSequencePercent(rows);
        drawSimpleViolin(xVel(i), y, velColors(i,:), violinWidth);
        if ~isempty(y)
            plot(xVel(i), median(y, 'omitnan'), 'o', 'MarkerFaceColor', velColors(i,:), ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');
        end
    end

    xline(numel(krValues) + 1, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, 'HandleVisibility', 'off');

    xticks([xKr xVel]);
    xticklabels([compose('%g', krValues) compose('%g', velocityValues)]);
    xlim([0.3, xVel(end) + 0.7]);
    ylim([0, 105]);
    ylabel('Buzz-readiness call (\% of full sequence)', 'Interpreter', 'latex');
    xlabel('$k_r \qquad\qquad\qquad\qquad\qquad v_0$', 'Interpreter', 'latex');
end

function plotDeltaTStarViolin(T, krValues, velocityValues, colorsKr, ceilingHitFrac, proxyConfig)
    hold on; box on; grid on; grid minor;

    xKr = 1:numel(krValues);
    xVel = (numel(krValues) + 2):(numel(krValues) + 1 + numel(velocityValues));
    violinWidth = 0.34;
    tbColor = [0.72 0.72 0.72];
    tbOffset = -0.16;
    dtOffset = 0.16;

    for i = 1:numel(krValues)
        rowsTb = abs(T.kr - krValues(i)) < 1e-12 & isfinite(T.(proxyConfig.VariableName));
        yTb = T.(proxyConfig.VariableName)(rowsTb);
        drawSimpleViolin(xKr(i) + tbOffset, yTb, tbColor, violinWidth * 0.80, 0.55, [0.25 0.25 0.25]);
        if ~isempty(yTb)
            plot(xKr(i) + tbOffset, median(yTb, 'omitnan'), 'o', 'MarkerFaceColor', 'k', ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 4, 'HandleVisibility', 'off');
        end

        rows = abs(T.kr - krValues(i)) < 1e-12 & isfinite(T.DeltaTStar_ms);
        y = T.DeltaTStar_ms(rows);
        drawSimpleViolin(xKr(i) + dtOffset, y, colorsKr(i,:), violinWidth * 0.68, 0.24);
        if ~isempty(y)
            plot(xKr(i) + dtOffset, median(y, 'omitnan'), 'o', 'MarkerFaceColor', colorsKr(i,:), ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');
        end
    end

    velColors = turbo(numel(velocityValues));
    for i = 1:numel(velocityValues)
        rows = abs(T.InitialVelocity_m_s - velocityValues(i)) < 1e-12 & isfinite(T.DeltaTStar_ms);
        y = T.DeltaTStar_ms(rows);
        drawSimpleViolin(xVel(i), y, velColors(i,:), violinWidth * 0.82, 0.24);
        if ~isempty(y)
            plot(xVel(i), median(y, 'omitnan'), 'o', 'MarkerFaceColor', velColors(i,:), ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');
        end
    end

    xline(numel(krValues) + 1, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, 'HandleVisibility', 'off');

    xticks([xKr xVel]);
    xticklabels([compose('%g', krValues) compose('%g', velocityValues)]);
    xlim([0.3, xVel(end) + 0.7]);
    ylabel(sprintf('$k_r$: $\\delta t^*$ at buzz readiness and %s (ms)', proxyConfig.AxisLabel), 'Interpreter', 'latex');
    yyaxis right
    ylim(yyaxisLimits(gca));
    yticks([]);
    ylabel('$v_0$: $\delta t^*$ at buzz readiness (ms)', 'Interpreter', 'latex');
    yyaxis left
    xlabel('$k_r \qquad\qquad\qquad\qquad\qquad v_0$', 'Interpreter', 'latex');

    % text(0.98, 0.97, sprintf('grey violins: %s', proxyConfig.LegendLabel), ...
    %     'Units', 'normalized', 'Interpreter', 'latex', ...
    %     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 6, ...
    %     'Color', [0.25 0.25 0.25]);
end

function plotReadinessTimeToContactViolin(T, krValues, velocityValues, colorsKr)
    hold on; box on; grid on; grid minor;

    xKr = 1:numel(krValues);
    xVel = (numel(krValues) + 2):(numel(krValues) + 1 + numel(velocityValues));
    violinWidth = 0.34;
    refColor = [0.72 0.72 0.72];
    leftOffset = -0.16;
    rightOffset = 0.16;

    for i = 1:numel(krValues)
        rowsRef = abs(T.kr - krValues(i)) < 1e-12 & isfinite(T.TimeFromDstopToTarget_ms);
        yRef = filterRobustOutliers(T.TimeFromDstopToTarget_ms(rowsRef));
        drawSimpleViolin(xKr(i) + leftOffset, yRef, refColor, violinWidth * 0.78, 0.55, [0.25 0.25 0.25]);
        if ~isempty(yRef)
            plot(xKr(i) + leftOffset, median(yRef, 'omitnan'), 'o', 'MarkerFaceColor', 'k', ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 4, 'HandleVisibility', 'off');
        end

        rows = abs(T.kr - krValues(i)) < 1e-12 & isfinite(T.TimeToTargetAtReadiness_ms);
        y = filterRobustOutliers(T.TimeToTargetAtReadiness_ms(rows));
        drawSimpleViolin(xKr(i) + rightOffset, y, colorsKr(i,:), violinWidth * 0.68, 0.24);
        if ~isempty(y)
            plot(xKr(i) + rightOffset, median(y, 'omitnan'), 'o', 'MarkerFaceColor', colorsKr(i,:), ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');
        end
    end

    velColors = turbo(numel(velocityValues));
    for i = 1:numel(velocityValues)
        rowsRef = abs(T.InitialVelocity_m_s - velocityValues(i)) < 1e-12 & isfinite(T.TimeFromDstopToTarget_ms);
        yRef = filterRobustOutliers(T.TimeFromDstopToTarget_ms(rowsRef));
        drawSimpleViolin(xVel(i) + leftOffset, yRef, refColor, violinWidth * 0.78, 0.55, [0.25 0.25 0.25]);
        if ~isempty(yRef)
            plot(xVel(i) + leftOffset, median(yRef, 'omitnan'), 'o', 'MarkerFaceColor', 'k', ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 4, 'HandleVisibility', 'off');
        end

        rows = abs(T.InitialVelocity_m_s - velocityValues(i)) < 1e-12 & isfinite(T.TimeToTargetAtReadiness_ms);
        y = filterRobustOutliers(T.TimeToTargetAtReadiness_ms(rows));
        drawSimpleViolin(xVel(i) + rightOffset, y, velColors(i,:), violinWidth * 0.68, 0.24);
        if ~isempty(y)
            plot(xVel(i) + rightOffset, median(y, 'omitnan'), 'o', 'MarkerFaceColor', velColors(i,:), ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');
        end
    end

    xline(numel(krValues) + 1, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, 'HandleVisibility', 'off');

    xticks([xKr xVel]);
    xticklabels([compose('%g', krValues) compose('%g', velocityValues)]);
    xlim([0.3, xVel(end) + 0.7]);
    set(gca, 'YScale', 'log');
    ylabel('Time to contact (ms)', 'Interpreter', 'latex');
    xlabel('$k_r \qquad\qquad\qquad\qquad\qquad v_0$', 'Interpreter', 'latex');
    text(0.98, 0.97, 'grey: $d_{\mathrm{stop}}\!\rightarrow\!$ target, colour: buzz readiness $\!\rightarrow\!$ target', ...
        'Units', 'normalized', 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 8, ...
        'Color', [0.25 0.25 0.25]);
end

function drawSimpleViolin(x0, y, faceColor, halfWidth, faceAlpha, edgeColor)
    if isempty(y)
        return
    end
    if nargin < 5
        faceAlpha = 0.18;
    end
    if nargin < 6
        edgeColor = faceColor;
    end

    y = y(isfinite(y));
    if numel(y) < 3
        scatter(repmat(x0, size(y)), y, 10, faceColor, 'filled', ...
            'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0.05, 'HandleVisibility', 'off');
        return
    end

    yGrid = linspace(min(y), max(y), 120);
    try
        f = ksdensity(y, yGrid, 'Function', 'pdf');
    catch
        scatter(repmat(x0, size(y)), y, 10, faceColor, 'filled', ...
            'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0.05, 'HandleVisibility', 'off');
        return
    end

    if all(~isfinite(f)) || max(f) <= 0
        return
    end

    w = halfWidth * (f / max(f));
    patch([x0 - w, fliplr(x0 + w)], [yGrid, fliplr(yGrid)], faceColor, ...
        'FaceAlpha', faceAlpha, 'EdgeColor', edgeColor, 'LineWidth', 1.15, 'HandleVisibility', 'off');
    plot([x0 x0], [prctile(y, 25), prctile(y, 75)], '-', 'Color', edgeColor, ...
        'LineWidth', 1.8, 'HandleVisibility', 'off');
end

function y = filterRobustOutliers(y)
    y = y(isfinite(y));
    if numel(y) < 5
        return
    end

    medY = median(y, 'omitnan');
    madY = mad(y, 1);
    if ~isfinite(madY) || madY == 0
        return
    end

    keep = abs(y - medY) <= 3 * 1.4826 * madY;
    y = y(keep);
end

function lims = yyaxisLimits(ax)
    lims = ax.YLim;
end

function makeExampleCurvePanel(curveTable, seqSummary, groupVar, groupValues, colors, keepMask, mainTitle, subTitle)
    nexttile;
    hold on; box on; grid on; grid minor;

    Skeep = seqSummary(keepMask & isfinite(seqSummary.BuzzReadinessDistance_m), :);
    for i = 1:numel(groupValues)
        gv = groupValues(i);
        rows = abs(Skeep.(groupVar) - gv) < 1e-12;
        Sg = Skeep(rows, :);
        if isempty(Sg)
            continue
        end

        [~, bestRow] = max(Sg.PrimeResponsivityNorm);
        seqID = Sg.SeqID(bestRow);
        C = curveTable(curveTable.SeqID == seqID & isfinite(curveTable.ResponsivityNorm), :);
        if isempty(C)
            continue
        end

        plot(C.ResponsivityDistance_m, C.ResponsivityNorm, '-', 'Color', colors(i,:), ...
            'LineWidth', 1.8, ...
            'DisplayName', makeGroupLabel(groupVar, gv));

        idxStar = unique(C.BuzzReadinessIndex);
        idxStar = idxStar(1);
        starRows = C.ResponsivityIndex == idxStar;
        plot(C.ResponsivityDistance_m(starRows), C.ResponsivityNorm(starRows), 'o', ...
            'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 6, ...
            'HandleVisibility', 'off');
    end

    yline(1, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xlabel('Anchor-distance midpoint (m)', 'Interpreter', 'latex');
    ylabel('Normalised responsivity, $\tilde{\mathcal{R}}_n$', 'Interpreter', 'latex');
    title({mainTitle, subTitle}, 'Interpreter', 'latex');
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 8, 'Box', 'on');
end

function label = makeGroupLabel(groupVar, value)
    switch char(groupVar)
        case 'kr'
            label = sprintf('$k_r=%g$', value);
        case 'InitialVelocity_m_s'
            label = sprintf('$v_0=%g$ m s$^{-1}$', value);
        case 'InitialCallDuration_ms'
            label = sprintf('$T_{c,0}=%g$ ms', value);
        otherwise
            label = sprintf('%s = %g', char(groupVar), value);
    end
end

function mask = fixedValueMask(T, varName, value)
    mask = abs(T.(varName) - value) < 1e-12;
end

function Tout = appendCompatible(Tin, Tnew)
    if isempty(Tin)
        Tout = Tnew;
        return
    end
    if isempty(Tnew)
        Tout = Tin;
        return
    end

    vTin = Tin.Properties.VariableNames;
    vTnew = Tnew.Properties.VariableNames;
    allVars = unique([vTin, vTnew], 'stable');

    for i = 1:numel(allVars)
        v = allVars{i};
        if ~ismember(v, vTin)
            Tin.(v) = missingColumn(height(Tin), Tnew.(v));
        end
        if ~ismember(v, vTnew)
            Tnew.(v) = missingColumn(height(Tnew), Tin.(v));
        end
    end

    Tnew = Tnew(:, Tin.Properties.VariableNames);
    Tout = [Tin; Tnew];
end

function x = missingColumn(n, template)
    if isnumeric(template)
        x = nan(n, 1);
    elseif islogical(template)
        x = false(n, 1);
    elseif isstring(template)
        x = strings(n, 1);
    else
        x = cell(n, 1);
    end
end
