%% Run multi-target SSG simulation and analysis from the responsivity core
% This script reproduces the earlier SSG-emergence workflow using the new
% simulateResponsivityCore framework. It:
%   1. simulates the original multi-target condition grid,
%   2. marks SSG patterns downstream from the raw timing table,
%   3. reproduces the compact event-level analyses and figures used in the
%      earlier analyse_ssg_patterns.m workflow.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir, 'fcn'));

resultsDir = fullfile(thisDir, 'results');
figureDir = fullfile(resultsDir, 'ssg_analysis_figures');

writeTables = false;
savePlots = false;  % switch on after figure tuning
runConfig = responsivityRunConfig();
if runConfig.OverrideOutputSwitches
    writeTables = runConfig.SaveTables;
    savePlots = runConfig.SaveFigures;
end
if runConfig.CloseFiguresBeforeRun
    close all;
end

if writeTables && ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
if savePlots && ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

%% Figure style
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');

%% Shared simulation grid (matched to the earlier SSG workflow)
initialDistances = 1:1:5;
meanVelocities = 2:2:10;
velocityJitterPc = 10:5:25;
nRepeats = 10;

baseParams = struct();
baseParams.kr = 5;
baseParams.maxCalls = 250;
baseParams.maxElapsedTime_s = 5;
baseParams.maxAnchorDistance_m = 10;
baseParams.interceptDistance_m = 0.15;
baseParams.initialCallDuration_s = 0.005;
baseParams.minCallDuration_s = 0.0005;

baseOpts = struct();
baseOpts.geometryMode = "3D";
baseOpts.initialTargetSeparation_m = 0.5;
baseOpts.timingMode = "motionAware";
baseOpts.phiMode = "constant";
baseOpts.phi = 0;
baseOpts.callDurationMode = "constant";
baseOpts.enforceMaxCallRate = false;
baseOpts.enforceResponseFloor = false;
baseOpts.batVelocityMode = "jittered";
baseOpts.batSteeringMode = "anchor";

conditions = [
    struct('name',"two_targets_stationary_anchor_nearest", 'numTargets',2, 'targetMotion',false, 'anchorMode',"nearest", 'code',"C2", 'displayOrder',2)
    struct('name',"two_targets_stationary_anchor_random",  'numTargets',2, 'targetMotion',false, 'anchorMode',"random",  'code',"C3", 'displayOrder',3)
    struct('name',"one_target_moving",                     'numTargets',1, 'targetMotion',true,  'anchorMode',"single",  'code',"C4", 'displayOrder',4)
    struct('name',"two_targets_moving_anchor_nearest",     'numTargets',2, 'targetMotion',true,  'anchorMode',"nearest", 'code',"C5", 'displayOrder',5)
    struct('name',"two_targets_moving_anchor_random",      'numTargets',2, 'targetMotion',true,  'anchorMode',"random",  'code',"C6", 'displayOrder',6)
    struct('name',"one_target_stationary",                 'numTargets',1, 'targetMotion',false, 'anchorMode',"single",  'code',"C1", 'displayOrder',1)
];

%% Downstream SSG definition
ssgParams = struct();
ssgParams.groupSizes = [2 3 4];
ssgParams.maxWithinCV = 0.05;
ssgParams.boundaryMultiplier = 1.20;

%% Simulate call tables
allCalls = table();
seqCounter = 0;
rng(1);

for cIdx = 1:numel(conditions)
    cond = conditions(cIdx);

    for d0 = initialDistances
    for vMean = meanVelocities
    for jitterPc = velocityJitterPc
    for rep = 1:nRepeats
        seqCounter = seqCounter + 1;

        params = baseParams;
        params.numSequences = 1;
        params.initialDistance_m = d0;
        params.initialBatSpeed_m_s = vMean;

        opts = baseOpts;
        opts.numTargets = cond.numTargets;
        opts.targetMotion = cond.targetMotion;
        opts.anchorMode = cond.anchorMode;
        opts.batVelocityJitterFrac = jitterPc / 100;

        if cond.targetMotion
            opts.targetVelocityMode = "stochastic";
            opts.targetVelocityScale = 0.5;
            opts.targetVelocityJitterFrac = (10 + 15 * rand()) / 100;
        else
            opts.targetVelocityMode = "stationary";
            opts.targetVelocityScale = 0;
            opts.targetVelocityJitterFrac = 0;
        end

        res = simulateResponsivityCore(params, opts);
        T = res.calls;
        if isempty(T)
            continue
        end

        T.SeqID(:) = seqCounter;
        T.ConditionName = repmat(cond.name, height(T), 1);
        T.ConditionID = repmat(cond.displayOrder, height(T), 1);
        T.ConditionCode = repmat(cond.code, height(T), 1);
        T.Repeat = repmat(rep, height(T), 1);
        T.InitialDistance_m = repmat(d0, height(T), 1);
        T.MeanVelocity_m_s = repmat(vMean, height(T), 1);
        T.BatVelocityJitter_percent = repmat(jitterPc, height(T), 1);
        T.IncludeTargetMotion = repmat(double(cond.targetMotion), height(T), 1);
        T.TwoTargetMode = repmat(double(cond.numTargets > 1), height(T), 1);
        T.AnchorSelectionMode = repmat(cond.anchorMode, height(T), 1);
        T.TargetVelocityScale = repmat(opts.targetVelocityScale, height(T), 1);
        T.TargetVelocityJitterFrac = repmat(opts.targetVelocityJitterFrac, height(T), 1);
        T = ensureTargetColumns(T, 2);

        if isempty(allCalls)
            allCalls = T;
        else
            T = alignToReferenceColumns(T, allCalls.Properties.VariableNames);
            allCalls = [allCalls; T]; %#ok<AGROW>
        end
    end
    end
    end
    end
end

%% Add SSG labels from the core timing output
[allCalls, eventTable, ssgSummary] = markSSGPatterns(allCalls, ssgParams);
allCalls = addSSGCompatibilityColumns(allCalls);
eventRows = buildSSGEventRows(allCalls);

%% Write tables
callFile = fullfile(resultsDir, 'core_simulated_temporal_patterns_ssg.csv');
eventFile = fullfile(resultsDir, 'core_ssg_event_summary.csv');
analysisEventFile = fullfile(resultsDir, 'core_ssg_event_analysis_summary.csv');

if writeTables
    writetable(allCalls, callFile);
    writetable(eventTable, eventFile);
    writetable(eventRows, analysisEventFile);
end

%% Basic setup
conditionCodeTable = table(string({conditions.code})', string({conditions.name})', ...
    [conditions.displayOrder]', ...
    'VariableNames', {'ConditionCode','ConditionName','DisplayOrder'});
conditionCodeTable = sortrows(conditionCodeTable, 'DisplayOrder');
conditionCodes = conditionCodeTable.ConditionCode;
conditionsStable = conditionCodeTable.ConditionName;
nConditions = height(conditionCodeTable);

if writeTables
    fprintf('Saved core SSG call table to: %s\n', callFile);
    fprintf('Saved core SSG event table to: %s\n', eventFile);
end
fprintf('Total calls: %d\n', height(allCalls));
fprintf('SSG calls: %d\n', sum(allCalls.IsSSG == 1));
fprintf('SSG events: %d\n', height(unique(allCalls(allCalls.IsSSG == 1, {'SeqID','SSG_ID'}))));

fprintf('\n=== Condition codes ===\n');
disp(conditionCodeTable);

conditionColours = [
    0.121 0.466 0.705
    1.000 0.498 0.054
    0.172 0.627 0.172
    0.839 0.153 0.157
    0.580 0.404 0.741
    0.549 0.337 0.294
    0.890 0.467 0.761
    0.498 0.498 0.498
];
if nConditions > size(conditionColours, 1)
    conditionColours = lines(nConditions);
end

%% Composite 2-by-2 SSG analysis figure
seqMap = table();
y = 0;
conditionBandStart = zeros(nConditions, 1);
conditionBandEnd = zeros(nConditions, 1);

for i = 1:nConditions
    cond = conditionsStable(i);
    seqIDs = unique(allCalls.SeqID(allCalls.ConditionName == cond));
    conditionBandStart(i) = y + 0.5;

    for s = 1:numel(seqIDs)
        y = y + 1;
        seqMap = [seqMap; table(cond, seqIDs(s), y, ...
            'VariableNames', {'ConditionName','SeqID','Y'})]; %#ok<AGROW>
    end

    conditionBandEnd(i) = y + 0.5;
end

totalSeq = y;
conditionBandMid = zeros(nConditions, 1);
for i = 1:nConditions
    conditionBandMid(i) = mean([conditionBandStart(i), conditionBandEnd(i)]);
end

isDoublet = eventRows.GroupSize == 2;
isTriplet = eventRows.GroupSize == 3;
isQuadruplet = eventRows.GroupSize == 4;
nDoublet = sum(isDoublet);
nTriplet = sum(isTriplet);
nQuadruplet = sum(isQuadruplet);
nSSGEvents = height(eventRows);
pDoublet = 100 * nDoublet / max(nSSGEvents, 1);
pTriplet = 100 * nTriplet / max(nSSGEvents, 1);
pQuadruplet = 100 * nQuadruplet / max(nSSGEvents, 1);

ssgConditions = unique(eventRows.ConditionName, 'stable');
histColours = lines(numel(ssgConditions));
binEdges = linspace(0, max(eventRows.AnchorDistance_m, [], 'omitnan'), 32);

conditionsWithSSG = unique(eventRows.ConditionName, 'stable');
sameCounts = zeros(numel(conditionsWithSSG), 1);
switchCounts = zeros(numel(conditionsWithSSG), 1);
for i = 1:numel(conditionsWithSSG)
    cond = conditionsWithSSG(i);
    sameCounts(i) = sum(eventRows.ConditionName == cond & eventRows.AnchorSameWithinSSG == true);
    switchCounts(i) = sum(eventRows.ConditionName == cond & eventRows.AnchorSwitchWithinSSG == true);
end

conditionCodesWithSSG = strings(numel(conditionsWithSSG), 1);
for i = 1:numel(conditionsWithSSG)
    conditionCodesWithSSG(i) = conditionCodeFromName(conditionsWithSSG(i), conditionsStable, conditionCodes);
end

figSSG = figure('Color', 'w', 'Position', [120 80 1000 940]);
tlSSG = tiledlayout(figSSG, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% A. Position of SSG events across sequences
axA = nexttile(tlSSG, 1);
hold(axA, 'on');

for i = 1:nConditions
    patch(axA, [0 100 100 0], ...
        [conditionBandStart(i) conditionBandStart(i) conditionBandEnd(i) conditionBandEnd(i)], ...
        conditionColours(i, :), ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 1.0);
end

for r = 1:height(seqMap)
    plot(axA, [0 100], [seqMap.Y(r) seqMap.Y(r)], '-', ...
        'Color', [0 0 0 0.08], 'LineWidth', 0.25);
end

for r = 1:height(eventRows)
    idxMap = seqMap.ConditionName == eventRows.ConditionName(r) & seqMap.SeqID == eventRows.SeqID(r);
    if any(idxMap)
        scatter(axA, eventRows.NormalisedPosition_percent(r), seqMap.Y(idxMap), ...
            16, 'filled', ...
            'MarkerFaceColor', [1.00 0.85 0.05], ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.25, ...
            'MarkerFaceAlpha', 1.0);
    end
end

hSSG = scatter(axA, nan, nan, 22, 'filled', ...
    'MarkerFaceColor', [1.00 0.85 0.05], 'MarkerEdgeColor', 'k');

xlabel(axA, 'Normalised sequence position (\%)', 'Interpreter', 'latex');
ylabel(axA, 'Simulation condition', 'Interpreter', 'latex');
title(axA, '\textbf{A. SSG position across simulated sequences}', 'Interpreter', 'latex');
xlim(axA, [0 100]);
ylim(axA, [0.5 totalSeq + 0.5]);
set(axA, 'YTick', conditionBandMid, 'YTickLabel', conditionCodes);
grid(axA, 'on'); grid(axA, 'minor'); box(axA, 'on');
pbaspect(axA, [1 1 1]);
lgdPosition = legend(axA, hSSG, {'Detected SSG'}, 'Location', 'southeast', 'Box', 'on', ...
    'FontSize', 8, 'Interpreter', 'latex');
hold(axA, 'off');

% B. Stable IPIs and their flanking intervals
axB = nexttile(tlSSG, 2);
hold(axB, 'on');
hStableDoublet = scatter(axB, eventRows.StableIPI_ms(isDoublet), eventRows.StableIPI_ms(isDoublet), ...
    24, 's', 'filled', ...
    'MarkerFaceColor', [0.05 0.60 0.20], ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.25, ...
    'MarkerFaceAlpha', 0.85);

hStableTriplet = scatter(axB, eventRows.StableIPI_ms(isTriplet), eventRows.StableIPI_ms(isTriplet), ...
    24, 's', 'filled', ...
    'MarkerFaceColor', [0.55 0.20 0.85], ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.25, ...
    'MarkerFaceAlpha', 0.85);

hStableQuadruplet = scatter(axB, eventRows.StableIPI_ms(isQuadruplet), eventRows.StableIPI_ms(isQuadruplet), ...
    28, 'd', 'filled', ...
    'MarkerFaceColor', [0.95 0.55 0.10], ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.25, ...
    'MarkerFaceAlpha', 0.85);

hLeft = scatter(axB, eventRows.StableIPI_ms, eventRows.LeftFlankIPI_ms, ...
    26, 'o', 'filled', ...
    'MarkerFaceColor', [0.05 0.30 0.95], ...
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', 0.62);

hRight = scatter(axB, eventRows.StableIPI_ms, eventRows.RightFlankIPI_ms, ...
    26, '^', 'filled', ...
    'MarkerFaceColor', [0.95 0.10 0.05], ...
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', 0.62);

maxVal = max([eventRows.StableIPI_ms; eventRows.LeftFlankIPI_ms; eventRows.RightFlankIPI_ms], [], 'omitnan');
hUnity = plot(axB, [0 maxVal], [0 maxVal], 'k--', 'LineWidth', 1.2);

xlabel(axB, 'Stable within-SSG IPI (ms)', 'Interpreter', 'latex');
ylabel(axB, 'Flanking IPI (ms)', 'Interpreter', 'latex');
title(axB, '\textbf{B. Stable SSG IPIs and flanking intervals}', 'Interpreter', 'latex');
lgdFlanks = legend(axB, [hStableDoublet hStableTriplet hStableQuadruplet hLeft hRight hUnity], ...
    {sprintf('Doublets: %d (%.1f\\%%)', nDoublet, pDoublet), ...
     sprintf('Triplets: %d (%.1f\\%%)', nTriplet, pTriplet), ...
     sprintf('Quadruplets: %d (%.1f\\%%)', nQuadruplet, pQuadruplet), ...
     'Left flank', ...
     'Right flank', ...
     'Unity line'}, ...
    'Location', 'southeast', 'Interpreter', 'latex', 'Box', 'on');
grid(axB, 'on'); grid(axB, 'minor'); box(axB, 'on');
pbaspect(axB, [1 1 1]);
hold(axB, 'off');

% C. Anchor-distance distributions
axC = nexttile(tlSSG, 3);
hold(axC, 'on');
for i = 1:numel(ssgConditions)
    cond = ssgConditions(i);
    idx = eventRows.ConditionName == cond;

    histogram(axC, eventRows.AnchorDistance_m(idx), binEdges, ...
        'Normalization', 'probability', ...
        'FaceColor', histColours(i, :), ...
        'EdgeColor', histColours(i, :), ...
        'FaceAlpha', 0.28, ...
        'LineWidth', 1.4, ...
        'DisplayName', conditionCodeFromName(cond, conditionsStable, conditionCodes));
end

xlabel(axC, 'Anchor target distance during SSG (m)', 'Interpreter', 'latex');
ylabel(axC, 'Probability', 'Interpreter', 'latex');
title(axC, '\textbf{C. Anchor distances at SSG occurrence}', 'Interpreter', 'latex');
lgdDistance = legend(axC, 'Location', 'northeast', 'Interpreter', 'latex', 'Box', 'on');
grid(axC, 'on'); grid(axC, 'minor'); box(axC, 'on');
pbaspect(axC, [1 1 1]);
hold(axC, 'off');

% D. Anchor identity within each SSG
axD = nexttile(tlSSG, 4);
b = bar(axD, categorical(conditionCodesWithSSG), [sameCounts switchCounts], 'grouped');
b(1).FaceColor = [0.20 0.70 0.30];
b(1).FaceAlpha = 0.35;
b(1).EdgeColor = [0.20 0.55 0.25];
b(2).FaceColor = [0.55 0.25 0.85];
b(2).FaceAlpha = 0.35;
b(2).EdgeColor = [0.45 0.15 0.75];

ylabel(axD, 'Number of SSG events', 'Interpreter', 'latex');
xlabel(axD, 'Simulation condition', 'Interpreter', 'latex');
title(axD, '\textbf{D. Within-group anchor consistency}', 'Interpreter', 'latex');
lgdConsistency = legend(axD, {'Same anchor within SSG', 'Anchor switch within SSG'}, ...
    'Location', 'northwest', 'Interpreter', 'latex');
grid(axD, 'on'); grid(axD, 'minor'); box(axD, 'on');
axD.XTickLabelRotation = 20;
pbaspect(axD, [1 1 1]);

formatLatex(figSSG, "full-square");
reserveLegendMetricMargin(lgdPosition, 0.08, axA);
reserveLegendMetricMargin(lgdFlanks, 0.10, axB);
reserveLegendMetricMargin(lgdDistance, 0.08, axC);
reserveLegendMetricMargin(lgdConsistency, 0.08, axD);

if savePlots
    exportPaperFigure(figSSG, fullfile(figureDir, 'ssg_analysis_composite'));
end

%% Print compact summaries
fprintf('\n=== SSG event summary by condition ===\n');

[G, condNames] = findgroups(eventRows.ConditionName);
nEvents = splitapply(@numel, eventRows.SSG_ID, G);
meanStableIPI = splitapply(@mean, eventRows.StableIPI_ms, G);
sdStableIPI = splitapply(@std, eventRows.StableIPI_ms, G);
modeStableIPI = splitapply(@mode, eventRows.StableIPI_ms, G);
q1StableIPI = splitapply(@(x) prctile(x, 25), eventRows.StableIPI_ms, G);
q3StableIPI = splitapply(@(x) prctile(x, 75), eventRows.StableIPI_ms, G);
meanAnchorDistance = splitapply(@mean, eventRows.AnchorDistance_m, G);
sdAnchorDistance = splitapply(@std, eventRows.AnchorDistance_m, G);
modeAnchorDistance = splitapply(@mode, eventRows.AnchorDistance_m, G);
q1AnchorDistance = splitapply(@(x) prctile(x, 25), eventRows.AnchorDistance_m, G);
q3AnchorDistance = splitapply(@(x) prctile(x, 75), eventRows.AnchorDistance_m, G);
sameAnchorPercent = 100 * splitapply(@mean, double(eventRows.AnchorSameWithinSSG), G);

condCodesForSummary = strings(numel(condNames), 1);
for i = 1:numel(condNames)
    condCodesForSummary(i) = conditionCodeFromName(condNames(i), conditionsStable, conditionCodes);
end

summaryByCondition = table(condCodesForSummary, condNames, nEvents, ...
    meanStableIPI, sdStableIPI, modeStableIPI, q1StableIPI, q3StableIPI, ...
    meanAnchorDistance, sdAnchorDistance, modeAnchorDistance, q1AnchorDistance, q3AnchorDistance, ...
    sameAnchorPercent, ...
    'VariableNames', {'ConditionCode','ConditionName','NumSSGEvents','MeanStableIPI_ms', ...
    'SDStableIPI_ms','ModeStableIPI_ms','Q1StableIPI_ms','Q3StableIPI_ms', ...
    'MeanAnchorDistance_m','SDAnchorDistance_m','ModeAnchorDistance_m','Q1AnchorDistance_m','Q3AnchorDistance_m', ...
    'PercentSameAnchorWithinSSG'});
summaryByCondition = sortrows(summaryByCondition, 'ConditionCode');

disp(summaryByCondition);

if writeTables
    summaryFile = fullfile(resultsDir, 'core_ssg_event_summary_by_condition.csv');
    writetable(summaryByCondition, summaryFile);
end

if savePlots
    fprintf('\nSaved figures in: %s\n', figureDir);
else
    fprintf('\nFigure saving is off: figures were displayed but not written to disk.\n');
end
if writeTables
    fprintf('Saved tables in: %s\n', resultsDir);
else
    fprintf('Table writing is off: analysis tables were not rewritten.\n');
end

%% Local helper functions
function T = alignToReferenceColumns(T, referenceNames)
    for i = 1:numel(referenceNames)
        name = referenceNames{i};
        if ~ismember(name, T.Properties.VariableNames)
            T.(name) = nan(height(T), 1);
        end
    end

    extraNames = setdiff(T.Properties.VariableNames, referenceNames, 'stable');
    orderedNames = [referenceNames(:); extraNames(:)]';
    T = T(:, orderedNames);
end

function T = addSSGCompatibilityColumns(T)
    T.ConditionName = string(T.ConditionName);
    T.AnchorSelectionMode = string(T.AnchorSelectionMode);
    T.StopReason = string(T.StopReason);
    T.IPI_ms = 1000 * T.IPI_s;
    T.AnchorDistanceToTarget_m = T.AnchorDistance_m;
    T.DistanceChange_m = T.DistanceChangeObserved_m;

    if ~ismember('Target2Distance_m', T.Properties.VariableNames)
        T.Target2Distance_m = nan(height(T), 1);
    end
end

function T = ensureTargetColumns(T, maxTargets)
    for k = 1:maxTargets
        names = {
            sprintf('Target%dDistance_m', k)
            sprintf('Target%dX_m', k)
            sprintf('Target%dY_m', k)
            sprintf('Target%dZ_m', k)
            sprintf('Target%dVx_m_s', k)
            sprintf('Target%dVy_m_s', k)
            sprintf('Target%dVz_m_s', k)
        };

        for ii = 1:numel(names)
            if ~ismember(names{ii}, T.Properties.VariableNames)
                T.(names{ii}) = nan(height(T), 1);
            end
        end
    end
end

function eventRows = buildSSGEventRows(T)
    ssgRows = T(T.IsSSG == 1, :);
    eventRows = table();

    if isempty(ssgRows)
        return
    end

    conditions = unique(T.ConditionName, 'stable');

    for i = 1:numel(conditions)
        cond = conditions(i);
        Tcond = T(T.ConditionName == cond, :);
        seqIDs = unique(Tcond.SeqID);

        for s = 1:numel(seqIDs)
            sid = seqIDs(s);
            S = Tcond(Tcond.SeqID == sid, :);
            nCalls = height(S);
            eventIDs = unique(S.SSG_ID(S.SSG_ID > 0));

            for e = 1:numel(eventIDs)
                eid = eventIDs(e);
                E = S(S.SSG_ID == eid, :);

                firstCall = min(E.CallNumber);
                lastCall = max(E.CallNumber);
                leftCall = firstCall - 1;
                rightCall = lastCall + 1;

                if leftCall < 1 || rightCall > nCalls
                    continue
                end

                leftRow = S(S.CallNumber == leftCall, :);
                rightRow = S(S.CallNumber == rightCall, :);

                newRow = table();
                newRow.ConditionName = cond;
                newRow.SeqID = sid;
                newRow.ConditionCode = E.ConditionCode(1);
                newRow.SSG_ID = eid;
                newRow.GroupSize = E.SSG_GroupSize_calls(1);
                newRow.FirstCall = firstCall;
                newRow.LastCall = lastCall;
                newRow.SequenceLength_calls = nCalls;
                newRow.NormalisedPosition_percent = 100 * mean([firstCall lastCall]) / nCalls;
                newRow.StableIPI_ms = E.SSG_WithinMedianIPI_ms(1);
                newRow.WithinCV = E.SSG_WithinCV(1);
                newRow.LeftFlankIPI_ms = leftRow.IPI_ms(1);
                newRow.RightFlankIPI_ms = rightRow.IPI_ms(1);
                newRow.LeftFlankRatio = E.SSG_LeftBoundaryRatio(1);
                newRow.RightFlankRatio = E.SSG_RightBoundaryRatio(1);
                newRow.AnchorDistance_m = mean(E.AnchorDistanceToTarget_m, 'omitnan');
                newRow.NearestDistance_m = mean(E.NearestTargetDistance_m, 'omitnan');

                anchors = unique(E.AnchorTargetID);
                newRow.AnchorSameWithinSSG = numel(anchors) == 1;
                newRow.AnchorSwitchWithinSSG = numel(anchors) > 1;

                eventRows = [eventRows; newRow]; %#ok<AGROW>
            end
        end
    end
end

function codeOut = conditionCodeFromName(conditionName, allConditionNames, allConditionCodes)
    idx = find(allConditionNames == string(conditionName), 1, 'first');
    if isempty(idx)
        codeOut = string(conditionName);
    else
        codeOut = allConditionCodes(idx);
    end
end
