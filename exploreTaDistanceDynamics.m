%% Explore acoustic acquisition interval versus distance
% Dynamic single-target simulations with motile targets and variable bat
% velocity. Three panels isolate the effects of k_r, phi, and approach
% speed on T_a as a function of effective anchor distance.
%
% Figures are displayed by default. Set savePlots = true after tuning.

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
showKrDistanceCallRateFigure = false;
showCrVrFigure = false;

if (savePlots || saveStats) && ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Shared simulation settings
paramsBase = struct();
paramsBase.c = 343;
paramsBase.kr = 5;
paramsBase.initialDistance_m = 10;
paramsBase.initialBatSpeed_m_s = 5;
paramsBase.initialCallDuration_s = 0.003;
paramsBase.minCallDuration_s = 0.0001;
paramsBase.maxCalls = 320;
paramsBase.maxElapsedTime_s = 10;
paramsBase.maxAnchorDistance_m = 10;
paramsBase.interceptDistance_m = 0.15;
paramsBase.numSequences = 1;

optsBase = struct();
optsBase.geometryMode = "3D";
optsBase.numTargets = 1;
optsBase.targetMotion = true;
optsBase.targetVelocityMode = "stochastic";
optsBase.targetVelocityScale = 0.5;
optsBase.targetVelocityJitterFrac = 0.12;
optsBase.batVelocityMode = "jittered";
optsBase.batVelocityJitterFrac = 0.12;
optsBase.anchorMode = "single";
optsBase.timingMode = "motionAware";
optsBase.callDurationMode = "previousTa";
optsBase.enforceMaxCallRate = true;
optsBase.callDurationJitter.enabled = true;
optsBase.callDurationJitter.mode = "additive";
optsBase.callDurationJitter.rho = 0.25;

nSequencesPerCondition = 75;
rng(33);

%% Condition sets
krValues = [3 5 7 9];
phiValues = [0 0.5 1];
velocityValues = [2 4 6 8 10];

fixedKr = 5;
fixedPhi = 0.5;
fixedVelocity_m_s = 5;

colorsKr = parula(numel(krValues));
colorsPhi = [
    0.00 0.45 0.70;  % blue
    0.90 0.62 0.00;  % orange
    0.00 0.62 0.45   % bluish green
];
colorsVelocity = turbo(numel(velocityValues));

%% Run simulations
Tkr = runConditionSweep(paramsBase, optsBase, "kr", krValues, ...
    fixedPhi, fixedVelocity_m_s, nSequencesPerCondition);
Tphi = runConditionSweep(paramsBase, optsBase, "phi", phiValues, ...
    fixedKr, fixedVelocity_m_s, nSequencesPerCondition);
Tvel = runConditionSweep(paramsBase, optsBase, "velocity", velocityValues, ...
    fixedKr, fixedPhi, nSequencesPerCondition);

fprintf('k_r sweep calls: %d\n', height(Tkr));
fprintf('phi sweep calls: %d\n', height(Tphi));
fprintf('velocity sweep calls: %d\n', height(Tvel));

%% Quantitative summaries
distanceBinEdges = linspace(0, 4.8, 19);

slopeSummary = table();
interactionSummary = table();
binnedSummary = table();

[slopeKrTa, interactionKrTa, binnedKrTa] = summariseSweep(Tkr, "Ta_ms", "kr", distanceBinEdges);
[slopeKrTb, interactionKrTb, binnedKrTb] = summariseSweep(Tkr, "Tb_ms", "kr", distanceBinEdges);
[slopePhiTa, interactionPhiTa, binnedPhiTa] = summariseSweep(Tphi, "Ta_ms", "phi", distanceBinEdges);
[slopePhiTb, interactionPhiTb, binnedPhiTb] = summariseSweep(Tphi, "Tb_ms", "phi", distanceBinEdges);
[slopeVelTa, interactionVelTa, binnedVelTa] = summariseSweep(Tvel, "Ta_ms", "velocity", distanceBinEdges);
[slopeVelTb, interactionVelTb, binnedVelTb] = summariseSweep(Tvel, "Tb_ms", "velocity", distanceBinEdges);

slopeSummary = [slopeKrTa; slopeKrTb; slopePhiTa; slopePhiTb; slopeVelTa; slopeVelTb];
interactionSummary = [interactionKrTa; interactionKrTb; interactionPhiTa; interactionPhiTb; interactionVelTa; interactionVelTb];
binnedSummary = [binnedKrTa; binnedKrTb; binnedPhiTa; binnedPhiTb; binnedVelTa; binnedVelTb];

phiConvergenceSummary = [ ...
    summarisePhiConvergence(Tphi, "Ta_ms", distanceBinEdges); ...
    summarisePhiConvergence(Tphi, "Tb_ms", distanceBinEdges)];

fprintf('\n=== Interaction tests ===\n');
disp(interactionSummary(:, {'SweepName','Response','NumRows','P_Distance','P_Condition','P_DistanceByCondition','R2'}));

fprintf('\n=== Slope summary (first rows) ===\n');
disp(slopeSummary(1:min(12,height(slopeSummary)), ...
    {'SweepName','Response','ConditionValue','NumRows','Slope','Intercept','RSquared'}));

fprintf('\n=== Phi convergence summary ===\n');
disp(phiConvergenceSummary(:, {'Response','PhiPair','NumBins','SlopeAbsDiffPerM','InterceptAbsDiff_ms','PValue','RSquared'}));

if saveStats
    writetable(slopeSummary, fullfile(outDir, 'TaTb_slope_summary.csv'));
    writetable(interactionSummary, fullfile(outDir, 'TaTb_interaction_tests.csv'));
    writetable(binnedSummary, fullfile(outDir, 'TaTb_binned_summary.csv'));
    writetable(phiConvergenceSummary, fullfile(outDir, 'TaTb_phi_convergence_summary.csv'));
end

%% Plot T_a and T_b versus distance
fig = figure('Color','w', 'Position', [80 80 1200 800]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing','compact', 'Padding','compact');
title(tl, '\textbf{Acoustic acquisition and behavioural intervals across distance in dynamic motile-target simulations}', ...
    'Interpreter','latex', 'FontSize', 15);

nexttile;
plotIntervalPanel(Tkr, "Ta_ms", krValues, colorsKr, ...
    @(v) sprintf('$k_r=%g$', v), ...
    '\textbf{A. $T_a$: responsivity coefficient}', ...
    sprintf('$\\phi=%.1f$, $v_0=%g$ m s$^{-1}$', fixedPhi, fixedVelocity_m_s), ...
    false, true);

nexttile;
plotIntervalPanel(Tphi, "Ta_ms", phiValues, colorsPhi, ...
    @(v) sprintf('$\\phi=%.1f$', v), ...
    '\textbf{B. $T_a$: echo-window choice}', ...
    sprintf('$k_r=%g$, $v_0=%g$ m s$^{-1}$', fixedKr, fixedVelocity_m_s), ...
    false, false);

nexttile;
plotIntervalPanel(Tvel, "Ta_ms", velocityValues, colorsVelocity, ...
    @(v) sprintf('$v_0=%g$ m s$^{-1}$', v), ...
    '\textbf{C. $T_a$: approach speed}', ...
    sprintf('$k_r=%g$, $\\phi=%.1f$', fixedKr, fixedPhi), ...
    false, false);

nexttile;
plotIntervalPanel(Tkr, "Tb_ms", krValues, colorsKr, ...
    @(v) sprintf('$k_r=%g$', v), ...
    '\textbf{D. $T_b$: responsivity coefficient}', ...
    sprintf('$\\phi=%.1f$, $v_0=%g$ m s$^{-1}$', fixedPhi, fixedVelocity_m_s), ...
    true, true);

nexttile;
plotIntervalPanel(Tphi, "Tb_ms", phiValues, colorsPhi, ...
    @(v) sprintf('$\\phi=%.1f$', v), ...
    '\textbf{E. $T_b$: echo-window choice}', ...
    sprintf('$k_r=%g$, $v_0=%g$ m s$^{-1}$', fixedKr, fixedVelocity_m_s), ...
    true, false);

nexttile;
plotIntervalPanel(Tvel, "Tb_ms", velocityValues, colorsVelocity, ...
    @(v) sprintf('$v_0=%g$ m s$^{-1}$', v), ...
    '\textbf{F. $T_b$: approach speed}', ...
    sprintf('$k_r=%g$, $\\phi=%.1f$', fixedKr, fixedPhi), ...
    true, false);

formatLatex(fig, "full-landscape");

if savePlots
    exportPaperFigure(fig, fullfile(outDir, 'explore_Ta_vs_distance_dynamic'));
end

%% Optional figure: k_r effect on distance-dependent call rate
if showKrDistanceCallRateFigure
    figKrCr = figure('Color','w', 'Position', [180 180 680 560]);
    hold on; box on; grid on; grid minor;
    for i = 1:numel(krValues)
        kr = krValues(i);
        rows = abs(Tkr.ConditionValue - kr) < 1e-12;
        S = Tkr(rows, :);

        scatter(S.AnchorDistance_m, S.CallRate_Hz, 8, colorsKr(i,:), ...
            'filled', ...
            'MarkerFaceAlpha', 0.75, ...
            'MarkerEdgeAlpha', 0.02, ...
            'DisplayName', sprintf('$k_r=%g$', kr));
    end

    xlim([0 4.8]);
    xlabel('Effective anchor distance, $d$ (m)', 'Interpreter','latex');
    ylabel('Call rate, $C_r$ (Hz)', 'Interpreter','latex');
    title({'\textbf{Responsivity coefficient shapes distance-dependent call rate}', ...
        sprintf('$\\phi=%.1f$, $v_0=%g$ m s$^{-1}$', fixedPhi, fixedVelocity_m_s)}, ...
        'Interpreter','latex');
    legend('Location','northeast', 'Interpreter','latex', 'FontSize', 9, 'Box','on');

    formatLatex(figKrCr, "half-square");

    if savePlots
        exportPaperFigure(figKrCr, fullfile(outDir, 'explore_Cr_vs_distance_by_kr_dynamic'));
    end
end

%% Optional separate figure: C_r versus radial velocity
if showCrVrFigure
    figCrVr = figure('Color','w', 'Position', [180 180 620 520]);
    plotCrVrPanel(Tvel, ...
        '\textbf{Call rate and radial velocity}', ...
        sprintf('$k_r=%g$, $\\phi=%.1f$', fixedKr, fixedPhi));

    formatLatex(figCrVr, "half-square");

    if savePlots
        exportPaperFigure(figCrVr, fullfile(outDir, 'explore_Cr_vs_vr_dynamic'));
    end
end

%% Local helper functions
function Tall = runConditionSweep(paramsBase, optsBase, sweepName, values, fixedA, fixedB, nSequences)
    Tall = table();
    seqCounter = 0;

    for i = 1:numel(values)
        value = values(i);

        for s = 1:nSequences
            seqCounter = seqCounter + 1;

            params = paramsBase;
            opts = optsBase;

            switch sweepName
                case "kr"
                    params.kr = value;
                    params.callDurationGain = 1 / params.kr;
                    opts.phi = fixedA;
                    params.initialBatSpeed_m_s = fixedB;
                    conditionValue = value;
                case "phi"
                    params.kr = fixedA;
                    params.callDurationGain = 1 / params.kr;
                    opts.phi = value;
                    params.initialBatSpeed_m_s = fixedB;
                    conditionValue = value;
                case "velocity"
                    params.kr = fixedA;
                    params.callDurationGain = 1 / params.kr;
                    opts.phi = fixedB;
                    params.initialBatSpeed_m_s = value;
                    conditionValue = value;
                otherwise
                    error('Unknown sweep name: %s', sweepName);
            end

            params.initialDistance_m = 5;
            params.maxCallRate_Hz = crMaxFromStop(params, params.initialBatSpeed_m_s);
            opts.rngSeed = [];

            res = simulateResponsivityCore(params, opts);
            T = res.calls;
            if isempty(T)
                continue
            end

            T.SeqID(:) = seqCounter;
            T.SweepName = repmat(string(sweepName), height(T), 1);
            T.ConditionValue = repmat(conditionValue, height(T), 1);
            T.kr = repmat(params.kr, height(T), 1);
            T.Phi = repmat(opts.phi, height(T), 1);
            T.InitialVelocity_m_s = repmat(params.initialBatSpeed_m_s, height(T), 1);
            T.Ta_ms = 1000 * T.Ta_s;
            T.Tdelay_ms = 1000 * T.Tdelay_s;
            T.Tphi_ms = 1000 * T.Tphi_s;
            T.Tb_ms = 1000 * T.TbEffective_s;

            Tall = appendCompatible(Tall, T);
        end
    end
end

function plotIntervalPanel(T, intervalName, values, colors, labelFcn, panelTitle, subtitleText, showXLabel, showYLabel)
    hold on; box on; grid on; grid minor;

    for i = 1:numel(values)
        value = values(i);
        rows = abs(T.ConditionValue - value) < 1e-12;
        S = T(rows, :);

        scatter(S.AnchorDistance_m, S.(intervalName), 8, colors(i,:), ...
            'filled', ...
            'MarkerFaceAlpha', 0.75, ...
            'MarkerEdgeAlpha', 0.02, ...
            'DisplayName', labelFcn(value));
    end

    xlim([0 4.8]);
    if showXLabel
        xlabel('Effective anchor distance, $d$ (m)', 'Interpreter','latex');
    else
        xlabel('');
    end
    if showYLabel
        if intervalName == "Ta_ms"
            ylabel('Acoustic acquisition interval, $T_a$ (ms)', 'Interpreter','latex');
        else
            ylabel('Behavioural interval, $T_b$ (ms)', 'Interpreter','latex');
        end
    else
        ylabel('');
    end
    title({panelTitle, subtitleText}, 'Interpreter','latex');
    legend('Location','northwest', 'Interpreter','latex', 'FontSize', 8, 'Box','on');
end

function plotCrVrPanel(T, panelTitle, subtitleText)
    hold on; box on; grid on; grid minor;

    validRows = isfinite(T.CallRate_Hz) & isfinite(T.RelativeClosingVelocity_m_s);
    S = T(validRows, :);
    keepRows = keepDbscanCore(S.CallRate_Hz, S.RelativeClosingVelocity_m_s);
    S = S(keepRows, :);

    scatter(S.CallRate_Hz, S.RelativeClosingVelocity_m_s, 11, [0.18 0.18 0.18], ...
        'filled', ...
        'MarkerFaceAlpha', 0.16, ...
        'MarkerEdgeAlpha', 0.02, ...
        'DisplayName','DBSCAN retained calls');

    xlabel('Call rate, $C_r$ (Hz)', 'Interpreter','latex');
    ylabel('Radial closing velocity, $v_r$ (m s$^{-1}$)', 'Interpreter','latex');
    title({panelTitle, subtitleText}, 'Interpreter','latex');
    legend('Location','northwest', 'Interpreter','latex', 'FontSize', 8, 'Box','on');
end

function [slopeTbl, interactionTbl, binnedTbl] = summariseSweep(T, responseVar, sweepName, edges)
    slopeTbl = fitPerConditionLines(T, responseVar, sweepName);
    interactionTbl = fitInteractionModel(T, responseVar, sweepName);
    binnedTbl = makeBinnedSummary(T, responseVar, sweepName, edges);
end

function slopeTbl = fitPerConditionLines(T, responseVar, sweepName)
    values = unique(T.ConditionValue);
    rowsOut = struct([]);
    for i = 1:numel(values)
        value = values(i);
        rows = abs(T.ConditionValue - value) < 1e-12 & ...
            isfinite(T.AnchorDistance_m) & isfinite(T.(responseVar));
        S = T(rows, :);
        if height(S) < 5
            continue
        end

        mdl = fitlm(S.AnchorDistance_m, S.(responseVar));
        coef = mdl.Coefficients;
        row = struct();
        row.SweepName = string(sweepName);
        row.Response = string(responseVar);
        row.ConditionValue = value;
        row.NumRows = height(S);
        row.Slope = coef.Estimate(2);
        row.SlopeCI_low = coef.Estimate(2) - 1.96 * coef.SE(2);
        row.SlopeCI_high = coef.Estimate(2) + 1.96 * coef.SE(2);
        row.Intercept = coef.Estimate(1);
        row.InterceptCI_low = coef.Estimate(1) - 1.96 * coef.SE(1);
        row.InterceptCI_high = coef.Estimate(1) + 1.96 * coef.SE(1);
        row.PSlope = coef.pValue(2);
        row.RSquared = mdl.Rsquared.Ordinary;
        if isempty(rowsOut)
            rowsOut = row;
        else
            rowsOut(end + 1, 1) = row; %#ok<AGROW>
        end
    end
    if isempty(rowsOut)
        slopeTbl = table();
    else
        slopeTbl = struct2table(rowsOut);
    end
end

function interactionTbl = fitInteractionModel(T, responseVar, sweepName)
    responseName = char(responseVar);
    rows = isfinite(T.AnchorDistance_m) & isfinite(T.(responseVar)) & isfinite(T.ConditionValue);
    S = T(rows, {'AnchorDistance_m', responseName, 'ConditionValue'});
    S.ResponseValue = S.(responseVar);
    S.ConditionGroup = categorical(string(S.ConditionValue));

    mdlDistance = fitlm(S, 'ResponseValue ~ AnchorDistance_m');
    mdlAdditive = fitlm(S, 'ResponseValue ~ AnchorDistance_m + ConditionGroup');
    mdlInteraction = fitlm(S, 'ResponseValue ~ AnchorDistance_m * ConditionGroup');

    pDistance = mdlDistance.Coefficients.pValue(strcmp(mdlDistance.CoefficientNames, 'AnchorDistance_m'));
    pCondition = nestedModelPValue(mdlDistance, mdlAdditive);
    pInteraction = nestedModelPValue(mdlAdditive, mdlInteraction);

    interactionTbl = table(string(sweepName), string(responseName), height(S), ...
        pDistance, pCondition, pInteraction, mdlInteraction.Rsquared.Ordinary, ...
        'VariableNames', {'SweepName','Response','NumRows','P_Distance','P_Condition','P_DistanceByCondition','R2'});
end

function pValue = nestedModelPValue(mdlReduced, mdlFull)
    sseReduced = mdlReduced.SSE;
    sseFull = mdlFull.SSE;
    dfeReduced = mdlReduced.DFE;
    dfeFull = mdlFull.DFE;

    dfNum = dfeReduced - dfeFull;
    dfDen = dfeFull;

    if dfNum <= 0 || dfDen <= 0 || sseReduced < sseFull
        pValue = NaN;
        return
    end

    msNum = (sseReduced - sseFull) / dfNum;
    msDen = sseFull / dfDen;

    if msDen <= 0
        pValue = NaN;
        return
    end

    fStat = msNum / msDen;
    pValue = 1 - fcdf(fStat, dfNum, dfDen);
end

function binnedTbl = makeBinnedSummary(T, responseVar, sweepName, edges)
    values = unique(T.ConditionValue);
    rowsOut = struct([]);
    binMid = edges(1:end-1) + diff(edges)/2;

    for i = 1:numel(values)
        value = values(i);
        useRows = abs(T.ConditionValue - value) < 1e-12 & ...
            isfinite(T.AnchorDistance_m) & isfinite(T.(responseVar));
        x = T.AnchorDistance_m(useRows);
        y = T.(responseVar)(useRows);

        for b = 1:numel(binMid)
            if b < numel(binMid)
                idx = x >= edges(b) & x < edges(b+1);
            else
                idx = x >= edges(b) & x <= edges(b+1);
            end
            if ~any(idx)
                continue
            end

            row = struct();
            row.SweepName = string(sweepName);
            row.Response = string(responseVar);
            row.ConditionValue = value;
            row.DistanceBinMid_m = binMid(b);
            row.N = nnz(idx);
            row.Mean_ms = mean(y(idx), 'omitnan');
            row.SD_ms = std(y(idx), 0, 'omitnan');
            row.Median_ms = median(y(idx), 'omitnan');
            row.Q1_ms = prctile(y(idx), 25);
            row.Q3_ms = prctile(y(idx), 75);

            if isempty(rowsOut)
                rowsOut = row;
            else
                rowsOut(end + 1, 1) = row; %#ok<AGROW>
            end
        end
    end

    if isempty(rowsOut)
        binnedTbl = table();
    else
        binnedTbl = struct2table(rowsOut);
    end
end

function phiTbl = summarisePhiConvergence(Tphi, responseVar, edges)
    phiPairs = [0 0.5; 0 1; 0.5 1];
    rowsOut = struct([]);
    binMid = edges(1:end-1) + diff(edges)/2;

    medByPhi = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for p = unique(Tphi.ConditionValue(:))'
        useRows = abs(Tphi.ConditionValue - p) < 1e-12 & ...
            isfinite(Tphi.AnchorDistance_m) & isfinite(Tphi.(responseVar));
        x = Tphi.AnchorDistance_m(useRows);
        y = Tphi.(responseVar)(useRows);
        med = nan(numel(binMid), 1);
        n = zeros(numel(binMid), 1);
        for b = 1:numel(binMid)
            if b < numel(binMid)
                idx = x >= edges(b) & x < edges(b+1);
            else
                idx = x >= edges(b) & x <= edges(b+1);
            end
            if any(idx)
                med(b) = median(y(idx), 'omitnan');
                n(b) = nnz(idx);
            end
        end
        medByPhi(num2str(p)) = struct('median', med, 'N', n);
    end

    for i = 1:size(phiPairs, 1)
        p1 = phiPairs(i, 1);
        p2 = phiPairs(i, 2);
        A = medByPhi(num2str(p1));
        B = medByPhi(num2str(p2));
        absDiff = abs(B.median - A.median);
        keep = isfinite(absDiff) & A.N > 0 & B.N > 0;
        if nnz(keep) < 3
            continue
        end

        mdl = fitlm(binMid(keep), absDiff(keep));
        coef = mdl.Coefficients;
        row = struct();
        row.Response = string(responseVar);
        row.PhiPair = sprintf('%.1f_vs_%.1f', p1, p2);
        row.NumBins = nnz(keep);
        row.SlopeAbsDiffPerM = coef.Estimate(2);
        row.InterceptAbsDiff_ms = coef.Estimate(1);
        row.PValue = coef.pValue(2);
        row.RSquared = mdl.Rsquared.Ordinary;

        if isempty(rowsOut)
            rowsOut = row;
        else
            rowsOut(end + 1, 1) = row; %#ok<AGROW>
        end
    end

    if isempty(rowsOut)
        phiTbl = table();
    else
        phiTbl = struct2table(rowsOut);
    end
end

function trend = fitSmoothTrend(x, y)
    nBins = 32;
    edges = linspace(min(x), max(x), nBins + 1);
    xBin = nan(nBins, 1);
    yBin = nan(nBins, 1);

    for i = 1:nBins
        rows = x >= edges(i) & x < edges(i + 1);
        if i == nBins
            rows = x >= edges(i) & x <= edges(i + 1);
        end
        if nnz(rows) >= 4
            xBin(i) = median(x(rows), 'omitnan');
            yBin(i) = median(y(rows), 'omitnan');
        end
    end

    keep = isfinite(xBin) & isfinite(yBin);
    xBin = xBin(keep);
    yBin = yBin(keep);

    if numel(xBin) >= 5
        ySmooth = smoothdata(yBin, 'movmedian', 5);
    else
        ySmooth = yBin;
    end

    trend = struct('x', xBin, 'y', ySmooth);
end

function keepRows = keepDbscanCore(x, y)
    keepRows = true(size(x));
    if numel(x) < 20
        return
    end

    X = [x(:), y(:)];
    X = (X - mean(X, 1, 'omitnan')) ./ std(X, 0, 1, 'omitnan');

    if exist('dbscan', 'file') ~= 2
        warning('DBSCAN is unavailable. Plotting all valid C_r--v_r points.');
        return
    end

    epsilon = 0.22;
    minPts = 18;
    labels = dbscan(X, epsilon, minPts);
    keepRows = labels ~= -1;

    if nnz(keepRows) < 0.25 * numel(keepRows)
        warning('DBSCAN removed too many points. Plotting all valid C_r--v_r points instead.');
        keepRows = true(size(x));
    end
end

function T = appendCompatible(T, newT)
    if isempty(T)
        T = newT;
        return
    end

    allNames = union(T.Properties.VariableNames, newT.Properties.VariableNames, 'stable');
    T = addMissingColumns(T, allNames);
    newT = addMissingColumns(newT, allNames);
    T = [T(:, allNames); newT(:, allNames)];
end

function T = addMissingColumns(T, allNames)
    for i = 1:numel(allNames)
        name = allNames{i};
        if ~ismember(name, T.Properties.VariableNames)
            T.(name) = nan(height(T), 1);
        end
    end
end

function crMax_Hz = crMaxFromStop(params, radialVelocity_m_s)
    crMax_Hz = (params.c + radialVelocity_m_s) / ...
        (2 * (1 + params.kr) * params.interceptDistance_m);
end
