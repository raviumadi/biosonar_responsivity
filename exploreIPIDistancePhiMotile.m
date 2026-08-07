%% Explore IPI versus distance under variable velocity and motile targets
% Motile-target, jittered-bat simulation with k_r = 5 and three phi
% conditions. The plot shows raw IPI and call rate from the simulated IPI.
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
showDurationDiagnostic = false;

if (savePlots || saveStats) && ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Simulation settings
phiValues = [0 0.5 1];
kr = 5;
nSequencesPerPhi = 100;

paramsBase = struct();
paramsBase.c = 343;
paramsBase.kr = kr;
paramsBase.initialDistance_m = 8;
paramsBase.initialBatSpeed_m_s = 5;
paramsBase.initialCallDuration_s = 0.005;
paramsBase.minCallDuration_s = 0.0001;
% Sequence safety cap; this is not the call-rate ceiling.
paramsBase.maxCalls = 220;
paramsBase.maxElapsedTime_s = 5;
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

rng(21);
allCalls = table();
seqCounter = 0;

for p = 1:numel(phiValues)
    phi = phiValues(p);

    for s = 1:nSequencesPerPhi
        seqCounter = seqCounter + 1;

        params = paramsBase;
        params.initialDistance_m = 5;
        params.initialBatSpeed_m_s = 4 + 5 * rand();
        params.maxCallRate_Hz = crMaxFromStop(params, params.initialBatSpeed_m_s);

        opts = optsBase;
        opts.rngSeed = [];
        opts.phi = phi;

        res = simulateResponsivityCore(params, opts);
        T = res.calls;
        if isempty(T)
            continue
        end

        T.SeqID(:) = seqCounter;
        T.Phi = repmat(phi, height(T), 1);
        T.InitialDistance_m = repmat(params.initialDistance_m, height(T), 1);
        T.MeanVelocity_m_s = repmat(params.initialBatSpeed_m_s, height(T), 1);
        T.MaxCallRate_Hz = repmat(params.maxCallRate_Hz, height(T), 1);
        T.IPI_ms = 1000 * T.IPI_s;
        T.Tcall_ms = 1000 * T.Tcall_s;
        T.CallRate_Hz = 1000 ./ T.IPI_ms;

        allCalls = appendCompatible(allCalls, T);
    end
end

fprintf('Simulated calls: %d\n', height(allCalls));
fprintf('Phi conditions: %s\n', strjoin(string(phiValues), ', '));
fprintf('Derived C_r,max range: %.2f--%.2f Hz\n', ...
    min(allCalls.MaxCallRate_Hz), max(allCalls.MaxCallRate_Hz));

%% Quantitative summaries for the phi-dependent timing patterns
distanceBinEdges = linspace(0, 4.8, 19);
bendWindow_m = [0.3 1.5];

ipiSlopeSummary = fitPerPhiLines(allCalls, "IPI_ms");
ipiInteractionSummary = fitPhiInteractionModel(allCalls, "IPI_ms");
ipiBinnedSummary = makePhiBinnedSummary(allCalls, "IPI_ms", distanceBinEdges);
ipiConvergenceSummary = summarisePhiConvergence(allCalls, "IPI_ms", distanceBinEdges);

crBinnedSummary = makePhiBinnedSummary(allCalls, "CallRate_Hz", distanceBinEdges);
crSpreadSummary = summariseCallRateSpread(allCalls, distanceBinEdges, bendWindow_m);
crSpreadInteractionSummary = fitSpreadInteractionModel(crSpreadSummary);

tcSlopeSummary = fitPerPhiLines(allCalls, "Tcall_ms");
tcInteractionSummary = fitPhiInteractionModel(allCalls, "Tcall_ms");
tcBinnedSummary = makePhiBinnedSummary(allCalls, "Tcall_ms", distanceBinEdges);
tcConvergenceSummary = summarisePhiConvergence(allCalls, "Tcall_ms", distanceBinEdges);

fprintf('\n=== Phi interaction tests ===\n');
disp([ipiInteractionSummary; tcInteractionSummary]);

fprintf('\n=== Call-rate spread summary in bend window %.2f--%.2f m ===\n', ...
    bendWindow_m(1), bendWindow_m(2));
disp(crSpreadInteractionSummary);

fprintf('\n=== IPI slope summary ===\n');
disp(ipiSlopeSummary);

fprintf('\n=== Call-duration slope summary ===\n');
disp(tcSlopeSummary);

fprintf('\n=== Phi convergence summary ===\n');
disp([ipiConvergenceSummary; tcConvergenceSummary]);

if saveStats
    writetable(ipiSlopeSummary, fullfile(outDir, 'phi_motile_IPI_slope_summary.csv'));
    writetable(ipiInteractionSummary, fullfile(outDir, 'phi_motile_IPI_interaction_tests.csv'));
    writetable(ipiBinnedSummary, fullfile(outDir, 'phi_motile_IPI_binned_summary.csv'));
    writetable(ipiConvergenceSummary, fullfile(outDir, 'phi_motile_IPI_convergence_summary.csv'));

    writetable(crBinnedSummary, fullfile(outDir, 'phi_motile_Cr_binned_summary.csv'));
    writetable(crSpreadSummary, fullfile(outDir, 'phi_motile_Cr_spread_summary.csv'));
    writetable(crSpreadInteractionSummary, fullfile(outDir, 'phi_motile_Cr_spread_tests.csv'));

    writetable(tcSlopeSummary, fullfile(outDir, 'phi_motile_Tc_slope_summary.csv'));
    writetable(tcInteractionSummary, fullfile(outDir, 'phi_motile_Tc_interaction_tests.csv'));
    writetable(tcBinnedSummary, fullfile(outDir, 'phi_motile_Tc_binned_summary.csv'));
    writetable(tcConvergenceSummary, fullfile(outDir, 'phi_motile_Tc_convergence_summary.csv'));
end

%% Plot IPI and call rate versus distance by phi
fig = figure('Color','w', 'Position', [120 160 1200 400]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');
title(tl, '\textbf{Timing consequences of echo-window choice under motile-target and variable-velocity conditions}, $k_r =5$', ...
    'Interpreter','latex', 'FontSize', 15);

phiColors = [
    0.00 0.45 0.70;  % blue
    0.90 0.62 0.00;  % orange
    0.00 0.62 0.45   % bluish green
];

nexttile;
hold on; box on; grid on; grid minor;
for p = 1:numel(phiValues)
    phi = phiValues(p);
    rows = allCalls.Phi == phi;
    S = allCalls(rows, :);

    scatter(S.AnchorDistance_m, S.IPI_ms, 8, phiColors(p,:), ...
        'filled', ...
        'MarkerFaceAlpha', 0.75, ...
        'MarkerEdgeAlpha', 0.04, ...
        'DisplayName', sprintf('$\\phi=%.1f$', phi));
end

xlim([0 4.8]);
xlabel('Effective anchor distance, $d$ (m)', 'Interpreter','latex');
ylabel('Interpulse interval, IPI (ms)', 'Interpreter','latex');
title('\textbf{A. IPI from sim. call timing}', 'Interpreter','latex');
legend('Location','northwest', 'Interpreter','latex', 'FontSize', 8, 'Box','on');

nexttile;
hold on; box on; grid on; grid minor;
for p = 1:numel(phiValues)
    phi = phiValues(p);
    rows = allCalls.Phi == phi;
    S = allCalls(rows, :);

    scatter(S.AnchorDistance_m, 1000 ./ S.IPI_ms, 8, phiColors(p,:), ...
        'filled', ...
        'MarkerFaceAlpha', 0.75, ...
        'MarkerEdgeAlpha', 0.04, ...
        'DisplayName', sprintf('$\\phi=%.1f$', phi));
end

xlim([0 4.8]);
xlabel('Effective anchor distance, $d$ (m)', 'Interpreter','latex');
ylabel('Call rate, $C_r$ (Hz)', 'Interpreter','latex');
title('\textbf{B. Call rate}', 'Interpreter','latex');
legend('Location','northeast', 'Interpreter','latex', 'FontSize', 8, 'Box','on');

nexttile;
hold on; box on; grid on; grid minor;
for p = 1:numel(phiValues)
    phi = phiValues(p);
    rows = allCalls.Phi == phi;
    S = allCalls(rows, :);

    scatter(S.AnchorDistance_m, 1000 * S.Tcall_s, 8, phiColors(p,:), ...
        'filled', ...
        'MarkerFaceAlpha', 0.75, ...
        'MarkerEdgeAlpha', 0.04, ...
        'DisplayName', sprintf('$\\phi=%.1f$', phi));

    fitRows = isfinite(S.AnchorDistance_m) & isfinite(S.Tcall_s);
    if nnz(fitRows) > 3
        xFit = S.AnchorDistance_m(fitRows);
        yFit = 1000 * S.Tcall_s(fitRows);
        fitObj = fitSmoothTrend(xFit, yFit);
        plot(fitObj.x, fitObj.y, '-', 'Color', phiColors(p,:), ...
            'LineWidth', 1.15, 'HandleVisibility','off');
    end
end

xlim([0 4.8]);
xlabel('Effective anchor distance, $d$ (m)', 'Interpreter','latex');
ylabel('Call duration, $T_c$ (ms)', 'Interpreter','latex');
title('\textbf{C. $T_c$ contraction with approach}', 'Interpreter','latex');
legend('Location','best', 'Interpreter','latex', 'FontSize', 8, 'Box','on');

formatLatex(fig, "full-wide");

if savePlots
    exportPaperFigure(fig, fullfile(outDir, 'explore_IPI_vs_distance_phi_motile'));
end

%% Optional diagnostic: does the 0.5 m feature mark call-duration contraction?
if showDurationDiagnostic
    dContract_m = (paramsBase.c + allCalls.MeanVelocity_m_s) * ...
        paramsBase.initialCallDuration_s / 2;
    fprintf('Call-duration contraction distance range: %.3f--%.3f m; median %.3f m\n', ...
        min(dContract_m), max(dContract_m), median(dContract_m, 'omitnan'));

    figDiag = figure('Color','w', 'Position', [190 190 1120 520]);
    tlDiag = tiledlayout(figDiag, 1, 2, 'TileSpacing','compact', 'Padding','compact');
    title(tlDiag, '\textbf{Diagnostic: call-duration contraction explains the near-field timing bend}', ...
        'Interpreter','latex', 'FontSize', 15);

    nexttile;
    hold on; box on; grid on; grid minor;
    for p = 1:numel(phiValues)
        phi = phiValues(p);
        rows = allCalls.Phi == phi;
        S = allCalls(rows, :);

        scatter(S.AnchorDistance_m, 1000 * S.Tcall_s, 12, phiColors(p,:), ...
            'filled', ...
            'MarkerFaceAlpha', 0.18, ...
            'MarkerEdgeAlpha', 0.04, ...
            'DisplayName', sprintf('$\\phi=%.1f$', phi));
    end
    xline(median(dContract_m, 'omitnan'), 'k:', '$d_{\\mathrm{contract}}$', ...
        'Interpreter','latex', 'LineWidth', 1.2, 'LabelOrientation','horizontal');
    xline(min(dContract_m), ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8, ...
        'HandleVisibility','off');
    xline(max(dContract_m), ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8, ...
        'HandleVisibility','off');
    yline(1000 * paramsBase.initialCallDuration_s, '--', '$T_{c,0}$', ...
        'Interpreter','latex', 'LineWidth', 1.0);
    yline(1000 * paramsBase.minCallDuration_s, '-.', '$T_{c,\min}$', ...
        'Interpreter','latex', 'LineWidth', 1.0);
    xlim([0 4.8]);
    xlabel('Effective anchor distance, $d$ (m)', 'Interpreter','latex');
    ylabel('Call duration, $T_c$ (ms)', 'Interpreter','latex');
    title('\textbf{A. Duration contraction begins near the bend}', 'Interpreter','latex');
    legend('Location','northeast', 'Interpreter','latex', 'FontSize', 8, 'Box','on');

    nexttile;
    hold on; box on; grid on; grid minor;
    for p = 1:numel(phiValues)
        phi = phiValues(p);
        rows = allCalls.Phi == phi;
        S = allCalls(rows, :);

        scatter(S.AnchorDistance_m, 1000 * S.Tphi_s, 12, phiColors(p,:), ...
            'filled', ...
            'MarkerFaceAlpha', 0.18, ...
            'MarkerEdgeAlpha', 0.04, ...
            'DisplayName', sprintf('$\\phi=%.1f$', phi));
    end
    xline(median(dContract_m, 'omitnan'), 'k:', '$d_{\\mathrm{contract}}$', ...
        'Interpreter','latex', 'LineWidth', 1.2, 'LabelOrientation','horizontal');
    xline(min(dContract_m), ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8, ...
        'HandleVisibility','off');
    xline(max(dContract_m), ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8, ...
        'HandleVisibility','off');
    xlim([0 4.8]);
    xlabel('Effective anchor distance, $d$ (m)', 'Interpreter','latex');
    ylabel('Selected echo-window contribution, $T_\phi$ (ms)', 'Interpreter','latex');
    title('\textbf{B. The $\phi$ contribution contracts with $T_c$}', 'Interpreter','latex');
    legend('Location','northeast', 'Interpreter','latex', 'FontSize', 8, 'Box','on');

    formatLatex(figDiag, "full-landscape");

    if savePlots
        exportPaperFigure(figDiag, fullfile(outDir, 'explore_call_duration_diagnostic_phi_motile'));
    end
end

%% Local helper functions
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

function trend = fitSmoothTrend(x, y)
    nBins = 28;
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

function slopeTbl = fitPerPhiLines(T, responseVar)
    phiLevels = unique(T.Phi(isfinite(T.Phi)));
    rowsOut = cell(numel(phiLevels), 1);

    for i = 1:numel(phiLevels)
        phi = phiLevels(i);
        rows = T.Phi == phi & isfinite(T.AnchorDistance_m) & isfinite(T.(responseVar));
        S = T(rows, {'AnchorDistance_m', char(responseVar)});
        S.Properties.VariableNames = {'Distance', 'Response'};

        mdl = fitlm(S, 'Response ~ Distance');
        rowsOut{i} = table(phi, height(S), mdl.Coefficients.Estimate(2), ...
            mdl.Coefficients.Estimate(1), mdl.Rsquared.Ordinary, ...
            'VariableNames', {'Phi', 'NumRows', 'Slope', 'Intercept', 'RSquared'});
    end

    slopeTbl = vertcat(rowsOut{:});
    slopeTbl.Response = repmat(string(responseVar), height(slopeTbl), 1);
    slopeTbl = movevars(slopeTbl, 'Response', 'Before', 'Phi');
end

function interactionTbl = fitPhiInteractionModel(T, responseVar)
    rows = isfinite(T.AnchorDistance_m) & isfinite(T.Phi) & isfinite(T.(responseVar));
    S = T(rows, {'AnchorDistance_m', char(responseVar), 'Phi'});
    S.Properties.VariableNames = {'Distance', 'Response', 'Phi'};
    S.PhiGroup = categorical(string(S.Phi));

    mdlDistance = fitlm(S, 'Response ~ Distance');
    mdlAdditive = fitlm(S, 'Response ~ Distance + PhiGroup');
    mdlInteraction = fitlm(S, 'Response ~ Distance * PhiGroup');

    pCondition = nestedModelPValue(mdlDistance, mdlAdditive);
    pInteraction = nestedModelPValue(mdlAdditive, mdlInteraction);
    pDistance = mdlInteraction.Coefficients.pValue(strcmp(mdlInteraction.CoefficientNames, 'Distance'));

    interactionTbl = table(string(responseVar), height(S), pDistance, pCondition, ...
        pInteraction, mdlInteraction.Rsquared.Ordinary, ...
        'VariableNames', {'Response', 'NumRows', 'P_Distance', ...
        'P_Phi', 'P_DistanceByPhi', 'R2'});
end

function binnedTbl = makePhiBinnedSummary(T, responseVar, edges)
    phiLevels = unique(T.Phi(isfinite(T.Phi)));
    rowsOut = {};

    for i = 1:numel(phiLevels)
        phi = phiLevels(i);
        phiRows = T.Phi == phi & isfinite(T.AnchorDistance_m) & isfinite(T.(responseVar));
        x = T.AnchorDistance_m(phiRows);
        y = T.(responseVar)(phiRows);

        for b = 1:numel(edges) - 1
            inBin = x >= edges(b) & x < edges(b + 1);
            if b == numel(edges) - 1
                inBin = x >= edges(b) & x <= edges(b + 1);
            end
            if nnz(inBin) < 5
                continue
            end

            yBin = y(inBin);
            q = quantile(yBin, [0.25 0.75]);
            p = quantile(yBin, [0.10 0.90]);
            rowsOut(end + 1, :) = { ...
                string(responseVar), phi, b, mean(edges(b:b+1)), ...
                edges(b), edges(b+1), nnz(inBin), ...
                mean(yBin, 'omitnan'), std(yBin, 'omitnan'), ...
                median(yBin, 'omitnan'), q(1), q(2), q(2) - q(1), ...
                p(1), p(2), p(2) - p(1)}; %#ok<AGROW>
        end
    end

    binnedTbl = cell2table(rowsOut, 'VariableNames', ...
        {'Response', 'Phi', 'BinIndex', 'DistanceBinMid_m', ...
        'DistanceBinStart_m', 'DistanceBinEnd_m', 'NumRows', ...
        'MeanValue', 'SDValue', 'MedianValue', 'Q1', 'Q3', 'IQR', ...
        'P10', 'P90', 'P10P90Width'});
end

function convergenceTbl = summarisePhiConvergence(T, responseVar, edges)
    binned = makePhiBinnedSummary(T, responseVar, edges);
    phiPairs = [0 0.5; 0 1; 0.5 1];
    rowsOut = {};

    for i = 1:size(phiPairs, 1)
        phiA = phiPairs(i, 1);
        phiB = phiPairs(i, 2);

        A = binned(binned.Phi == phiA, {'DistanceBinMid_m', 'MedianValue'});
        B = binned(binned.Phi == phiB, {'DistanceBinMid_m', 'MedianValue'});
        AB = innerjoin(A, B, 'Keys', 'DistanceBinMid_m');
        AB.AbsDiff = abs(AB.MedianValue_A - AB.MedianValue_B);

        mdl = fitlm(AB, 'AbsDiff ~ DistanceBinMid_m');
        rowsOut(end + 1, :) = { ...
            string(responseVar), sprintf('%.1f_vs_%.1f', phiA, phiB), ...
            height(AB), abs(mdl.Coefficients.Estimate(2)), ...
            abs(mdl.Coefficients.Estimate(1)), ...
            mdl.Coefficients.pValue(2), mdl.Rsquared.Ordinary}; %#ok<AGROW>
    end

    convergenceTbl = cell2table(rowsOut, 'VariableNames', ...
        {'Response', 'PhiPair', 'NumBins', 'SlopeAbsDiffPerM', ...
        'InterceptAbsDiff', 'PValue', 'RSquared'});
end

function spreadTbl = summariseCallRateSpread(T, edges, bendWindow)
    binned = makePhiBinnedSummary(T, "CallRate_Hz", edges);
    binned.IsBendWindow = binned.DistanceBinMid_m >= bendWindow(1) & ...
        binned.DistanceBinMid_m <= bendWindow(2);
    spreadTbl = binned;
end

function spreadSummary = fitSpreadInteractionModel(spreadTbl)
    bendTbl = spreadTbl(spreadTbl.IsBendWindow, :);
    bendTbl.PhiGroup = categorical(string(bendTbl.Phi));

    mdlDistance = fitlm(bendTbl, 'IQR ~ DistanceBinMid_m');
    mdlAdditive = fitlm(bendTbl, 'IQR ~ DistanceBinMid_m + PhiGroup');
    mdlInteraction = fitlm(bendTbl, 'IQR ~ DistanceBinMid_m * PhiGroup');

    phiLevels = unique(bendTbl.Phi);
    rowsOut = cell(numel(phiLevels), 1);
    for i = 1:numel(phiLevels)
        phi = phiLevels(i);
        S = bendTbl(bendTbl.Phi == phi, :);
        [peakIQR, idx] = max(S.IQR);
        rowsOut{i} = table(phi, height(S), peakIQR, S.DistanceBinMid_m(idx), ...
            mean(S.IQR, 'omitnan'), mean(S.P10P90Width, 'omitnan'), ...
            nestedModelPValue(mdlDistance, mdlAdditive), ...
            nestedModelPValue(mdlAdditive, mdlInteraction), ...
            mdlInteraction.Rsquared.Ordinary, ...
            'VariableNames', {'Phi', 'NumBinsInBendWindow', 'PeakIQR_Hz', ...
            'DistanceAtPeakIQR_m', 'MeanIQR_Bend_Hz', 'MeanP10P90Width_Bend_Hz', ...
            'P_Phi', 'P_DistanceByPhi', 'R2'});
    end

    spreadSummary = vertcat(rowsOut{:});
end

function pValue = nestedModelPValue(mdlReduced, mdlFull)
    sseReduced = mdlReduced.SSE;
    sseFull = mdlFull.SSE;
    dfeReduced = mdlReduced.DFE;
    dfeFull = mdlFull.DFE;

    dfNum = dfeReduced - dfeFull;
    dfDen = dfeFull;

    if dfNum <= 0 || dfDen <= 0 || sseFull <= 0
        pValue = NaN;
        return
    end

    fStat = ((sseReduced - sseFull) / dfNum) / (sseFull / dfDen);
    if fStat < 0
        pValue = 1;
    else
        pValue = 1 - fcdf(fStat, dfNum, dfDen);
    end
end
