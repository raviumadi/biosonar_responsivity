%% Explore processed field-call data against the responsivity framework
% First-pass exploratory script for the processed field dataset.
% The goal here is to import the CSV cleanly, inspect the basic timing and
% kinematic variables, and flag suspect rows using robust sequence-wise
% criteria before deriving framework-specific quantities downstream.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'fcn'));
dataFile = fullfile(thisDir, 'data', 'vof_processed_data.csv');

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
outDir = fullfile(thisDir, 'results', 'field_data_exploration');
showFlaggedPoints = false;

% Reproducibility cache. When true, compatible saved simulation tables are
% loaded instead of rerunning the expensive stochastic ensembles. Set this
% false to force regeneration; newly generated tables are still saved when
% saveSimulationData is true.
loadSavedSimulationData = true;
saveSimulationData = true;
if runConfig.OverrideOutputSwitches
    saveSimulationData = runConfig.SaveTables;
end
simulationCacheDir = fullfile(outDir, 'simulation_cache');

% Manuscript figures (enabled by default)
makeObservableComparisonFigure = true;
makeTimingDiagnosticsFigure = true;

% Optional exploratory figures (disabled by default)
makeDurationEnvelopeTestFigure = false;
makeLegacyTimingDecompositionFigure = false;
makeIndependentDiagnosticFigure = false;
makeExtendedPhiProfileFigure = false;
makeD0PhiProfileFigure = false;
makeDistanceBinnedEchoWindowFigure = false;

% The d0-by-phi sensitivity supports numerical results in the manuscript,
% but its heatmap is optional and controlled separately above.
runD0PhiSensitivity = true;

matchSimulationVelocityToField = true;
velocityMatchBinWidth_m_s = 0.5;
velocityMatchSeed = 11;
diagnosticMinVelocityBinCount = 20;
profileVelocityMatchBinWidth_m_s = 1;
profileVelocityLimits_m_s = [2 9];
matchDurationDistanceSampleCountToField = true;
mainSimulationSeedBase = 42000;
phiProfileValues = 0:0.05:2;
phiOriginalSimulationLimit = 1;
phiProfileNumSequencesPerPhi = 100;
phiProfileQuantiles = 0.01:0.01:0.99;
phiProfileSimulationSeedBase = 84000;
phiDistanceBinWidth_m = 0.5;
d0ProfileValues_m = 2:1:8;
d0PhiProfileValues = 0:0.1:2;
d0PhiNumSequencesPerCombination = 30;
d0PhiSimulationSeedBase = 84000;
defaultSimulationInitialDistance_m = 5;
krModel = 5;
cSound = 343;
maxVelocityForKeep_m_s = 10;
velocityColourLimits_m_s = [0 maxVelocityForKeep_m_s];
minDistanceForHyperbolaFit_m = 0.05;
kinematicVelocityCandidates_m_s = [2 4 6 8 10];
simPhiValues = [0 0.5 1];
simNumSequencesPerPhi = 100;
durationCompatibilityNumReplicates = 20;
durationCompatibilitySeedBase = 51000;
durationCompatibilitySeedStep = 1000;
minCallRateForDurationFit_Hz = 5;
applyFieldDurationCapForDurationFit = true;
maxFieldDurationForDurationFit_ms = 10;
durationEnvelopeNumBins = 18;
durationEnvelopeQuantile = 0.10;
durationMedianQuantile = 0.50;

madThresh = 5;
jumpMadThresh = 5;

rng(42);

if (savePlots || saveSimulationData) && ~exist(outDir, 'dir')
    mkdir(outDir);
end
if saveSimulationData && ~exist(simulationCacheDir, 'dir')
    mkdir(simulationCacheDir);
end

%% Figure style
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');

% Shared field--simulation colour convention used throughout the primary
% comparison figures: warm orange for field data and cyan-blue for simulations.
fieldColour = [0.90 0.50 0.05];
simColour = [0.00 0.55 0.85];

%% Load and standardise
T = readtable(dataFile);
T = sortrows(T, {'SeqNum', 'Timestamp'});

% Derived convenience variables
T.CallRate_Hz = T.Rate;
T.IPI_s = T.IPI;
T.Velocity_m_s = T.Velocity;
T.Distance_m = T.Distance;
T.Duration_s = T.Duration;
T.Level_dB = T.Level;
T.DurationDiff_s = T.DurationDiff;
T.Timestamp_s = T.Timestamp;

T.IPI_ms = 1000 * T.IPI_s;
T.Duration_ms = 1000 * T.Duration_s;
T.DurationDiff_ms = 1000 * T.DurationDiff_s;
T.krModel = repmat(krModel, height(T), 1);
T.Ta_fromIPI_s = T.IPI_s ./ (1 + krModel);
T.Tb_fromIPI_s = (krModel ./ (1 + krModel)) .* T.IPI_s;
T.Ta_fromIPI_ms = 1000 * T.Ta_fromIPI_s;
T.Tb_fromIPI_ms = 1000 * T.Tb_fromIPI_s;

%% Flag physically implausible values
T.IsPhysicalOutlier = false(height(T), 1);
T.IsPhysicalOutlier = T.IsPhysicalOutlier | isnan(T.CallRate_Hz) | T.CallRate_Hz <= 0;
T.IsPhysicalOutlier = T.IsPhysicalOutlier | isnan(T.IPI_s) | T.IPI_s <= 0;
T.IsPhysicalOutlier = T.IsPhysicalOutlier | isnan(T.Velocity_m_s) | T.Velocity_m_s < 0;
T.IsPhysicalOutlier = T.IsPhysicalOutlier | isnan(T.Duration_s) | T.Duration_s <= 0;
T.IsPhysicalOutlier = T.IsPhysicalOutlier | isnan(T.Timestamp_s);

%% Flag robust sequence-wise outliers
T.IsRateOutlier = false(height(T), 1);
T.IsIPIOutlier = false(height(T), 1);
T.IsVelocityOutlier = false(height(T), 1);
T.IsDurationOutlier = false(height(T), 1);
T.IsDurationDiffOutlier = false(height(T), 1);

T.IsRateJumpOutlier = false(height(T), 1);
T.IsIPIJumpOutlier = false(height(T), 1);
T.IsVelocityJumpOutlier = false(height(T), 1);
T.IsDurationJumpOutlier = false(height(T), 1);

seqIDs = unique(T.SeqNum);

for i = 1:numel(seqIDs)
    idx = T.SeqNum == seqIDs(i);
    S = T(idx, :);

    T.IsRateOutlier(idx) = robustFlag(S.CallRate_Hz, madThresh);
    T.IsIPIOutlier(idx) = robustFlag(S.IPI_ms, madThresh);
    T.IsVelocityOutlier(idx) = robustFlag(S.Velocity_m_s, madThresh);
    T.IsDurationOutlier(idx) = robustFlag(S.Duration_ms, madThresh);
    T.IsDurationDiffOutlier(idx) = robustFlag(S.DurationDiff_ms, madThresh);

    T.IsRateJumpOutlier(idx) = robustJumpFlag(S.CallRate_Hz, jumpMadThresh);
    T.IsIPIJumpOutlier(idx) = robustJumpFlag(S.IPI_ms, jumpMadThresh);
    T.IsVelocityJumpOutlier(idx) = robustJumpFlag(S.Velocity_m_s, jumpMadThresh);
    T.IsDurationJumpOutlier(idx) = robustJumpFlag(S.Duration_ms, jumpMadThresh);

end

T.IsAnyOutlier = T.IsPhysicalOutlier | T.IsRateOutlier | T.IsIPIOutlier | ...
    T.IsVelocityOutlier | T.IsDurationOutlier | T.IsDurationDiffOutlier | ...
    T.IsRateJumpOutlier | T.IsIPIJumpOutlier | T.IsVelocityJumpOutlier | ...
    T.IsDurationJumpOutlier;

T.IsKeep = ~T.IsAnyOutlier;
T.IsKeep = T.IsKeep & T.Velocity_m_s <= maxVelocityForKeep_m_s;

%% Print basic summary
fprintf('Loaded field dataset: %s\n', dataFile);
fprintf('Rows: %d\n', height(T));
fprintf('Sequences: %d\n', numel(unique(T.SeqNum)));
fprintf('Time span: %.3f--%.3f s\n', min(T.Timestamp_s), max(T.Timestamp_s));
fprintf('Call-rate range: %.2f--%.2f Hz\n', min(T.CallRate_Hz), max(T.CallRate_Hz));
fprintf('Velocity range: %.2f--%.2f m/s\n', min(T.Velocity_m_s), max(T.Velocity_m_s));
fprintf('Distance range: %.3f--%.3f m\n', min(T.Distance_m), max(T.Distance_m));
fprintf('Duration range: %.3f--%.3f ms\n', min(T.Duration_ms), max(T.Duration_ms));
fprintf('Flagged rows: %d (%.2f%%)\n', sum(T.IsAnyOutlier), 100 * mean(T.IsAnyOutlier));
fprintf('Kept-row velocity cap: <= %.2f m/s\n', maxVelocityForKeep_m_s);
fprintf('Timing decomposition uses fixed k_r = %.2f\n', krModel);
if applyFieldDurationCapForDurationFit
    fprintf('Duration-fit cap applied to field data: <= %.2f ms\n', maxFieldDurationForDurationFit_ms);
else
    fprintf('Duration-fit cap applied to field data: none\n');
end

fitIdx = T.IsKeep & isfinite(T.Distance_m) & isfinite(T.CallRate_Hz) & ...
    T.Distance_m >= minDistanceForHyperbolaFit_m & T.CallRate_Hz > 0;
fprintf('\nKinematic call-rate reference curves against inter-call displacement\n');
fprintf('(Using C_r = v / Delta x with fixed velocity references)\n');
disp(table(kinematicVelocityCandidates_m_s(:), 'VariableNames', {'velocity_m_s'}));

%% Build a simulation dataset and match its call-level velocity distribution
fieldVelocityForMatching_m_s = T.Velocity_m_s(T.IsKeep & ...
    isfinite(T.Velocity_m_s) & T.Velocity_m_s >= 0);
simulationInitialSpeedRange_m_s = [ ...
    max(0.1, min(fieldVelocityForMatching_m_s, [], 'omitnan')), ...
    maxVelocityForKeep_m_s];
mainSimulationConfig = struct( ...
    'CacheVersion', 1, ...
    'Kr', krModel, ...
    'MaxVelocity_m_s', maxVelocityForKeep_m_s, ...
    'PhiValues', simPhiValues, ...
    'NumSequencesPerPhi', simNumSequencesPerPhi, ...
    'InitialSpeedRange_m_s', simulationInitialSpeedRange_m_s, ...
    'SeedBase', mainSimulationSeedBase, ...
    'InitialDistance_m', defaultSimulationInitialDistance_m);
mainSimulationCacheFile = fullfile(simulationCacheDir, ...
    'field_validation_main.mat');
[S, mainSimulationLoaded] = loadSimulationCache( ...
    mainSimulationCacheFile, mainSimulationConfig, ...
    loadSavedSimulationData);
if ~mainSimulationLoaded
    S = runFieldValidationSimulation(thisDir, krModel, ...
        maxVelocityForKeep_m_s, simPhiValues, simNumSequencesPerPhi, ...
        simulationInitialSpeedRange_m_s, mainSimulationSeedBase, ...
        defaultSimulationInitialDistance_m);
end
mainSimulationNeedsCacheWrite = ~mainSimulationLoaded || ...
    width(S) > numel(fieldValidationCacheVariables());
S = trimFieldValidationSimulationTable(S);
if mainSimulationNeedsCacheWrite
    saveSimulationCache(mainSimulationCacheFile, S, ...
        mainSimulationConfig, saveSimulationData);
end
velocityMatchEdges_m_s = 0:velocityMatchBinWidth_m_s:maxVelocityForKeep_m_s;
if velocityMatchEdges_m_s(end) < maxVelocityForKeep_m_s
    velocityMatchEdges_m_s(end + 1) = maxVelocityForKeep_m_s;
end
if matchSimulationVelocityToField
    [S, velocityMatchSummary] = matchSimulationVelocityDistribution( ...
        S, fieldVelocityForMatching_m_s, velocityMatchEdges_m_s, velocityMatchSeed);
else
    velocityMatchSummary = table();
end
fprintf('\nValidation simulation\n');
fprintf('Main simulation source: %s\n', ...
    simulationSourceLabel(mainSimulationLoaded));
fprintf('Simulated calls kept: %d\n', height(S));
fprintf('Simulation phi values: %s\n', strjoin(string(simPhiValues), ', '));
fprintf('Simulation velocity-matching switch: %d\n', matchSimulationVelocityToField);
if matchSimulationVelocityToField
    fprintf('Velocity-matching bins: %.2f m/s\n', velocityMatchBinWidth_m_s);
    disp(velocityMatchSummary);
end

%% Quantitative comparison for rate-displacement structure
fieldRateRows = T.IsKeep & isfinite(T.CallRate_Hz) & isfinite(T.Distance_m) & ...
    T.CallRate_Hz > 0 & T.Distance_m > 0;
simRateRows = isfinite(S.CallRate_Hz) & isfinite(S.InterCallDisplacement_m) & ...
    S.CallRate_Hz > 0 & S.InterCallDisplacement_m > 0;

fieldEffectiveVelocity_m_s = T.CallRate_Hz(fieldRateRows) .* T.Distance_m(fieldRateRows);
simEffectiveVelocity_m_s = S.CallRate_Hz(simRateRows) .* S.InterCallDisplacement_m(simRateRows);
fieldEffectiveDistance_m = T.Distance_m(fieldRateRows);
simEffectiveDistance_m = S.InterCallDisplacement_m(simRateRows);

fieldBandIdx = assignNearestBand(fieldEffectiveVelocity_m_s, kinematicVelocityCandidates_m_s);
simBandIdx = assignNearestBand(simEffectiveVelocity_m_s, kinematicVelocityCandidates_m_s);

fieldBandCounts = accumarray(fieldBandIdx, 1, [numel(kinematicVelocityCandidates_m_s), 1], @sum, 0);
simBandCounts = accumarray(simBandIdx, 1, [numel(kinematicVelocityCandidates_m_s), 1], @sum, 0);

fieldBandProps = fieldBandCounts / sum(fieldBandCounts);
simBandProps = simBandCounts / sum(simBandCounts);

[chi2Stat, chi2Df, chi2PValue] = chiSquareHomogeneity([fieldBandCounts.'; simBandCounts.']);

velocityBandSummary = table( ...
    kinematicVelocityCandidates_m_s(:), ...
    fieldBandCounts(:), fieldBandProps(:), ...
    simBandCounts(:), simBandProps(:), ...
    'VariableNames', {'VelocityBand_m_s','FieldCount','FieldProportion','SimulationCount','SimulationProportion'});

velocityOverallSummary = table( ...
    ["Field"; "Simulation"], ...
    [numel(fieldEffectiveVelocity_m_s); numel(simEffectiveVelocity_m_s)], ...
    [median(fieldEffectiveVelocity_m_s, 'omitnan'); median(simEffectiveVelocity_m_s, 'omitnan')], ...
    [prctile(fieldEffectiveVelocity_m_s, 25); prctile(simEffectiveVelocity_m_s, 25)], ...
    [prctile(fieldEffectiveVelocity_m_s, 75); prctile(simEffectiveVelocity_m_s, 75)], ...
    [mean(fieldEffectiveVelocity_m_s, 'omitnan'); mean(simEffectiveVelocity_m_s, 'omitnan')], ...
    'VariableNames', {'Dataset','NumRows','MedianEffectiveVelocity_m_s', ...
    'Q1EffectiveVelocity_m_s','Q3EffectiveVelocity_m_s','MeanEffectiveVelocity_m_s'});

fprintf('\nEffective-velocity comparison from call-rate and inter-call displacement\n');
disp(velocityOverallSummary);
fprintf('Velocity-band occupancy comparison: chi-square = %.3f, df = %d, p = %.6g\n', ...
    chi2Stat, chi2Df, chi2PValue);
disp(velocityBandSummary);

vEffBinEdges = linspace(min([fieldEffectiveDistance_m; simEffectiveDistance_m], [], 'omitnan'), ...
    max([fieldEffectiveDistance_m; simEffectiveDistance_m], [], 'omitnan'), 16);
[vEffFieldCtr, vEffFieldQ25, vEffFieldQ50, vEffFieldQ75, vEffFieldN] = ...
    binnedQuantiles(fieldEffectiveDistance_m, fieldEffectiveVelocity_m_s, vEffBinEdges, [25 50 75]);
[vEffSimCtr, vEffSimQ25, vEffSimQ50, vEffSimQ75, vEffSimN] = ...
    binnedQuantiles(simEffectiveDistance_m, simEffectiveVelocity_m_s, vEffBinEdges, [25 50 75]);

vEffFieldSummary = table(vEffFieldCtr(:), vEffFieldQ25(:), vEffFieldQ50(:), vEffFieldQ75(:), vEffFieldN(:), ...
    'VariableNames', {'DistanceMid_m','FieldQ25_m_s','FieldMedian_m_s','FieldQ75_m_s','FieldNumRows'});
vEffSimSummary = table(vEffSimCtr(:), vEffSimQ25(:), vEffSimQ50(:), vEffSimQ75(:), vEffSimN(:), ...
    'VariableNames', {'DistanceMid_m','SimulationQ25_m_s','SimulationMedian_m_s','SimulationQ75_m_s','SimulationNumRows'});

%% Timing-structure comparison: field vs simulation
timingComparison = buildTimingComparisonTable(T, S);
timingInteractionSummary = table();
timingResidualSummary = table();
timingResidualBinSummary = table();

responseSpecs = { ...
    'IPI_ms', 'IPI_ms'; ...
    'Ta_fromIPI_ms', 'Ta_ms'; ...
    'Tb_fromIPI_ms', 'Tb_ms'};

for iResp = 1:size(responseSpecs, 1)
    fieldVar = responseSpecs{iResp, 1};
    simVar = responseSpecs{iResp, 2};
    responseLabel = erase(fieldVar, '_fromIPI');

    summaryRow = fitFieldSimulationInteractionModel(timingComparison, fieldVar, simVar, responseLabel);
    timingInteractionSummary = appendCompatible(timingInteractionSummary, summaryRow);

    residualSummaryRows = summariseResidualSpread(timingComparison, fieldVar, simVar, responseLabel);
    timingResidualSummary = appendCompatible(timingResidualSummary, residualSummaryRows);

    residualBinRows = summariseResidualSpreadByDistance(timingComparison, fieldVar, simVar, responseLabel);
    timingResidualBinSummary = appendCompatible(timingResidualBinSummary, residualBinRows);
end

fprintf('\nField-simulation interaction tests for displacement- and velocity-dependent timing\n');
disp(timingInteractionSummary);
fprintf('\nResidual spread summary after removing distance and velocity structure\n');
disp(timingResidualSummary);

fieldDurationFitIdx = T.IsKeep;
if applyFieldDurationCapForDurationFit
    fieldDurationFitIdx = fieldDurationFitIdx & T.Duration_ms <= maxFieldDurationForDurationFit_ms;
end

fieldDurationFit = fitDurationRateHyperbola(T.CallRate_Hz(fieldDurationFitIdx), ...
    T.Duration_ms(fieldDurationFitIdx), minCallRateForDurationFit_Hz);
simDurationFit = fitDurationRateHyperbola(S.CallRate_Hz, S.Tcall_ms, minCallRateForDurationFit_Hz);

fprintf('\nDuration-vs-rate hyperbolic fit\n');
disp(table( ...
    ["Field"; "Simulation"], ...
    [fieldDurationFit.NumRows; simDurationFit.NumRows], ...
    [fieldDurationFit.Offset_ms; simDurationFit.Offset_ms], ...
    [fieldDurationFit.Scale_msHz; simDurationFit.Scale_msHz], ...
    [fieldDurationFit.RMSE_ms; simDurationFit.RMSE_ms], ...
    [fieldDurationFit.RSquared; simDurationFit.RSquared], ...
    'VariableNames', {'Dataset','NumRows','Offset_ms','Scale_msHz','RMSE_ms','RSquared'}));
if applyFieldDurationCapForDurationFit
    fprintf('Field rows excluded from duration fit by duration cap: %d\n', ...
        sum(T.IsKeep & T.Duration_ms > maxFieldDurationForDurationFit_ms));
end

fieldDurationEnvelopeFit = fitDurationRateEnvelope( ...
    T.CallRate_Hz(fieldDurationFitIdx), T.Duration_ms(fieldDurationFitIdx), ...
    minCallRateForDurationFit_Hz, durationEnvelopeNumBins, durationEnvelopeQuantile);
simDurationEnvelopeFit = fitDurationRateEnvelope( ...
    S.CallRate_Hz, S.Tcall_ms, minCallRateForDurationFit_Hz, ...
    durationEnvelopeNumBins, durationEnvelopeQuantile);
fieldDurationMedianFit = fitDurationRateEnvelope( ...
    T.CallRate_Hz(fieldDurationFitIdx), T.Duration_ms(fieldDurationFitIdx), ...
    minCallRateForDurationFit_Hz, durationEnvelopeNumBins, durationMedianQuantile);
simDurationMedianFit = fitDurationRateEnvelope( ...
    S.CallRate_Hz, S.Tcall_ms, minCallRateForDurationFit_Hz, ...
    durationEnvelopeNumBins, durationMedianQuantile);

%% Duration-distance compatibility envelope for field durations
fieldEnvelopeRows = T.IsKeep & isfinite(T.Ta_fromIPI_ms) & T.Ta_fromIPI_ms > 0 & ...
    isfinite(T.Duration_ms) & T.Duration_ms > 0 & ...
    isfinite(T.CallRate_Hz) & isfinite(T.Distance_m);
Fd = T(fieldEnvelopeRows, :);
vrProxy_m_s = Fd.Velocity_m_s;
rangeFactor = 0.5 * (cSound + vrProxy_m_s);
Fd.DistanceFar_m = rangeFactor .* Fd.Ta_fromIPI_s;
Fd.DistanceNear_m = rangeFactor .* max(Fd.Ta_fromIPI_s - Fd.Duration_s, 0);
Fd.DistanceMid_m = 0.5 * (Fd.DistanceNear_m + Fd.DistanceFar_m);
Fd.DistanceEnvelopeWidth_m = Fd.DistanceFar_m - Fd.DistanceNear_m;

distKeep = isfinite(Fd.DistanceNear_m) & isfinite(Fd.DistanceFar_m) & ...
    Fd.DistanceFar_m >= Fd.DistanceNear_m & isfinite(Fd.Duration_ms);
Fd = Fd(distKeep, :);

durationCompatibilitySeeds = durationCompatibilitySeedBase + ...
    durationCompatibilitySeedStep .* ...
    (0:(durationCompatibilityNumReplicates - 1));
durationSimulationConfig = struct( ...
    'CacheVersion', 2, ...
    'Kr', krModel, ...
    'MaxVelocity_m_s', maxVelocityForKeep_m_s, ...
    'PhiValues', simPhiValues, ...
    'NumSequencesPerPhi', simNumSequencesPerPhi, ...
    'Seeds', durationCompatibilitySeeds);
durationSimulationCacheFile = fullfile(simulationCacheDir, ...
    'duration_distance_replicates.mat');
[durationSimulationReplicates, durationSimulationLoaded] = ...
    loadSimulationCache(durationSimulationCacheFile, ...
    durationSimulationConfig, loadSavedSimulationData);
if ~durationSimulationLoaded
    durationSimulationReplicates = cell( ...
        durationCompatibilityNumReplicates, 1);
    fprintf('\nGenerating %d seeded duration--distance ensembles\n', ...
        durationCompatibilityNumReplicates);
    for iReplicate = 1:durationCompatibilityNumReplicates
        durationSimulationReplicates{iReplicate} = ...
            runDurationDistanceSimulation(thisDir, krModel, ...
            maxVelocityForKeep_m_s, simPhiValues, ...
            simNumSequencesPerPhi, ...
            durationCompatibilitySeeds(iReplicate));
    end
end
durationSimulationNeedsCacheWrite = ~durationSimulationLoaded;
for iReplicate = 1:durationCompatibilityNumReplicates
    durationSimulationNeedsCacheWrite = ...
        durationSimulationNeedsCacheWrite || ...
        width(durationSimulationReplicates{iReplicate}) > ...
        numel(durationDistanceCacheVariables());
    durationSimulationReplicates{iReplicate} = ...
        trimDurationDistanceSimulationTable( ...
        durationSimulationReplicates{iReplicate});
end
if durationSimulationNeedsCacheWrite
    saveSimulationCache(durationSimulationCacheFile, ...
        durationSimulationReplicates, durationSimulationConfig, ...
        saveSimulationData);
end

compatPctReplicates = nan(durationCompatibilityNumReplicates, 1);
durationEnvelopeResults = cell(durationCompatibilityNumReplicates, 1);
for iReplicate = 1:durationCompatibilityNumReplicates
    durationEnvelopeResults{iReplicate} = ...
        evaluateDurationEnvelopeCompatibility(Fd, ...
        durationSimulationReplicates{iReplicate}, ...
        matchDurationDistanceSampleCountToField, ...
        durationCompatibilitySeeds(iReplicate) + 17, 35);
    compatPctReplicates(iReplicate) = ...
        durationEnvelopeResults{iReplicate}.CompatibilityPct;
end
compatPctMedian = median(compatPctReplicates, 'omitnan');
compatPctQ1 = prctile(compatPctReplicates, 25);
compatPctQ3 = prctile(compatPctReplicates, 75);
[~, representativeCompatibilityIdx] = min( ...
    abs(compatPctReplicates - compatPctMedian));
representativeEnvelope = ...
    durationEnvelopeResults{representativeCompatibilityIdx};
Sdist = representativeEnvelope.SimulationTable;
distBinCtr = representativeEnvelope.DistanceBinCentre_m;
simQ05 = representativeEnvelope.SimulationQ05_ms;
simQ50 = representativeEnvelope.SimulationMedian_ms;
simQ95 = representativeEnvelope.SimulationQ95_ms;
simN = representativeEnvelope.SimulationNumCalls;
Fd.IsEnvelopeCompatible = ...
    representativeEnvelope.IsFieldCompatible;

durationCompatibilitySummary = table( ...
    durationCompatibilityNumReplicates, compatPctMedian, ...
    compatPctQ1, compatPctQ3, ...
    durationCompatibilitySeeds(representativeCompatibilityIdx), ...
    compatPctReplicates(representativeCompatibilityIdx), ...
    'VariableNames', {'NumSeededReplicates', ...
    'MedianCompatible_pct', 'Q1Compatible_pct', 'Q3Compatible_pct', ...
    'RepresentativeSeed', 'RepresentativeCompatible_pct'});
fprintf('\nDuration--distance envelope compatibility across seeded ensembles\n');
fprintf('Simulation source: %s\n', ...
    simulationSourceLabel(durationSimulationLoaded));
disp(durationCompatibilitySummary);

fprintf('\nDuration-vs-rate lower-envelope hyperbolic fit\n');
disp(table( ...
    ["Field"; "Simulation"], ...
    [fieldDurationEnvelopeFit.NumRows; simDurationEnvelopeFit.NumRows], ...
    [fieldDurationEnvelopeFit.NumBinsUsed; simDurationEnvelopeFit.NumBinsUsed], ...
    [fieldDurationEnvelopeFit.Offset_ms; simDurationEnvelopeFit.Offset_ms], ...
    [fieldDurationEnvelopeFit.Scale_msHz; simDurationEnvelopeFit.Scale_msHz], ...
    [fieldDurationEnvelopeFit.RMSE_ms; simDurationEnvelopeFit.RMSE_ms], ...
    [fieldDurationEnvelopeFit.RSquared; simDurationEnvelopeFit.RSquared], ...
    'VariableNames', {'Dataset','NumRows','NumBinsUsed','Offset_ms','Scale_msHz','RMSE_ms','RSquared'}));

fprintf('\nDuration-vs-rate binwise-median hyperbolic fit\n');
disp(table( ...
    ["Field"; "Simulation"], ...
    [fieldDurationMedianFit.NumRows; simDurationMedianFit.NumRows], ...
    [fieldDurationMedianFit.NumBinsUsed; simDurationMedianFit.NumBinsUsed], ...
    [fieldDurationMedianFit.Offset_ms; simDurationMedianFit.Offset_ms], ...
    [fieldDurationMedianFit.Scale_msHz; simDurationMedianFit.Scale_msHz], ...
    [fieldDurationMedianFit.RMSE_ms; simDurationMedianFit.RMSE_ms], ...
    [fieldDurationMedianFit.RSquared; simDurationMedianFit.RSquared], ...
    'VariableNames', {'Dataset','NumRows','NumBinsUsed','Offset_ms','Scale_msHz','RMSE_ms','RSquared'}));

%% Field versus simulation: primary observable comparison
fieldRateXLim = [min(T.Distance_m(T.IsKeep), [], 'omitnan') max(T.Distance_m(T.IsKeep), [], 'omitnan')];
fieldRateYLim = [min(T.CallRate_Hz(T.IsKeep), [], 'omitnan') max(T.CallRate_Hz(T.IsKeep), [], 'omitnan')];
fieldIpiXLim = [min(T.Distance_m(T.IsKeep), [], 'omitnan') max(T.Distance_m(T.IsKeep), [], 'omitnan')];
fieldIpiYLim = [min(T.IPI_ms(T.IsKeep), [], 'omitnan') max(T.IPI_ms(T.IsKeep), [], 'omitnan')];
fieldDurXLim = [min(T.CallRate_Hz(T.IsKeep), [], 'omitnan') max(T.CallRate_Hz(T.IsKeep), [], 'omitnan')];
fieldDurYLim = [min(T.Duration_ms(T.IsKeep), [], 'omitnan') max(T.Duration_ms(T.IsKeep), [], 'omitnan')];

if makeObservableComparisonFigure
fig1 = figure('Color', 'w', 'Position', [100 100 1000 600]);
tl1 = tiledlayout(fig1, 2, 6, 'Padding', 'compact', 'TileSpacing', 'compact');

axA = nexttile(tl1, [1 2]);
curveColours = lines(numel(kinematicVelocityCandidates_m_s));
distRateHandles = gobjects(numel(kinematicVelocityCandidates_m_s), 1);
distRateLabels = strings(numel(kinematicVelocityCandidates_m_s), 1);
pointVelocityIdx = nan(sum(T.IsKeep), 1);
keepVelocities = T.Velocity_m_s(T.IsKeep);
for ii = 1:numel(keepVelocities)
    [~, pointVelocityIdx(ii)] = min(abs(kinematicVelocityCandidates_m_s - keepVelocities(ii)));
end
hold(axA, 'on');
for ii = 1:numel(kinematicVelocityCandidates_m_s)
    idxLocal = false(height(T), 1);
    tmp = find(T.IsKeep);
    idxLocal(tmp(pointVelocityIdx == ii)) = true;
    scatter(axA, T.Distance_m(idxLocal), T.CallRate_Hz(idxLocal), 10, curveColours(ii, :), 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75);
end
if showFlaggedPoints
    scatter(axA, T.Distance_m(T.IsAnyOutlier), T.CallRate_Hz(T.IsAnyOutlier), 14, ...
        [0.85 0.10 0.10], 'x', 'LineWidth', 0.8);
end
xlabel(axA, 'Inter-call displacement (m)', 'Interpreter', 'latex');
ylabel(axA, 'Call rate (Hz)', 'Interpreter', 'latex');
title(axA, '\textbf{A. Field call rate}', 'Interpreter', 'latex');
xFit = linspace(min(T.Distance_m(fitIdx)), max(T.Distance_m(fitIdx)), 300);
for ii = 1:numel(kinematicVelocityCandidates_m_s)
    vRef = kinematicVelocityCandidates_m_s(ii);
    yFit = vRef ./ xFit;
    distRateHandles(ii) = plot(axA, xFit, yFit, '-', 'Color', curveColours(ii, :), 'LineWidth', 1.8);
    distRateLabels(ii) = sprintf('$v = %.0f$ m s$^{-1}$', vRef);
end
legend(axA, distRateHandles, distRateLabels, 'Location', 'northeast');
hold(axA, 'off');
xlim(axA, fieldRateXLim);
ylim(axA, fieldRateYLim);
grid(axA, 'on'); grid(axA, 'minor'); box(axA, 'on');

axB = nexttile(tl1, [1 2]);
plotDistanceRatePanel(axB, S, kinematicVelocityCandidates_m_s, curveColours, ...
    minDistanceForHyperbolaFit_m, '\textbf{B. Simulated call rate}');
xlim(axB, fieldRateXLim);
ylim(axB, fieldRateYLim);

axC = nexttile(tl1, [1 2]);
hold(axC, 'on');
fieldBand = isfinite(vEffFieldCtr) & isfinite(vEffFieldQ25) & isfinite(vEffFieldQ50) & isfinite(vEffFieldQ75) & vEffFieldN > 0;
simBand = isfinite(vEffSimCtr) & isfinite(vEffSimQ25) & isfinite(vEffSimQ50) & isfinite(vEffSimQ75) & vEffSimN > 0;
if any(fieldBand)
    patch(axC, ...
        [vEffFieldCtr(fieldBand); flipud(vEffFieldCtr(fieldBand))], ...
        [vEffFieldQ25(fieldBand); flipud(vEffFieldQ75(fieldBand))], ...
        [0.92 0.74 0.14], 'FaceAlpha', 0.28, 'EdgeColor', 'none');
    plot(axC, vEffFieldCtr(fieldBand), vEffFieldQ50(fieldBand), '-', ...
        'Color', fieldColour, 'LineWidth', 2.0);
end
if any(simBand)
    patch(axC, ...
        [vEffSimCtr(simBand); flipud(vEffSimCtr(simBand))], ...
        [vEffSimQ25(simBand); flipud(vEffSimQ75(simBand))], ...
        [0.20 0.60 0.90], 'FaceAlpha', 0.22, 'EdgeColor', 'none');
    plot(axC, vEffSimCtr(simBand), vEffSimQ50(simBand), '--', ...
        'Color', simColour, 'LineWidth', 2.0);
end
hold(axC, 'off');
xlabel(axC, 'Inter-call displacement (m)', 'Interpreter', 'latex');
ylabel(axC, '$v_{\mathrm{eff}} = C_r \Delta x$ (m s$^{-1}$)', 'Interpreter', 'latex');
title(axC, '\textbf{C. Effective velocity}', ...
    'Interpreter', 'latex');
legend(axC, {'field IQR', 'field median', 'sim. IQR', 'sim. median'}, ...
    'Location', 'best');
xlim(axC, fieldRateXLim);
grid(axC, 'on'); grid(axC, 'minor'); box(axC, 'on');

axD = nexttile(tl1, [1 3]);
hold(axD, 'on');
scatter(axD, S.CallRate_Hz, S.Tcall_ms, 10, [0.08 0.18 0.45], 'filled', ...
    'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0.35);
scatter(axD, T.CallRate_Hz(T.IsKeep), T.Duration_ms(T.IsKeep), 10, ...
    [0.90 0.72 0.12], 'filled', 'MarkerFaceAlpha', 0.55, 'MarkerEdgeAlpha', 0.55);
plotDurationFitCurve(axD, fieldDurationMedianFit, fieldDurXLim, '--', 1.8, [0.90 0.40 0.05]);
plotDurationFitCurve(axD, simDurationMedianFit, fieldDurXLim, '--', 1.8, [0.00 0.75 0.85]);
hold(axD, 'off');
xlabel(axD, 'Call rate (Hz)', 'Interpreter', 'latex');
ylabel(axD, 'Call duration (ms)', 'Interpreter', 'latex');
title(axD, '\textbf{D. Call duration--rate relation}', 'Interpreter', 'latex');
xlim(axD, fieldDurXLim);
ylim(axD, fieldDurYLim);
lgdDuration = legend(axD, {'field', 'simulation', 'field binwise-median fit', 'sim. binwise-median fit'}, ...
    'Location', 'northeast', 'Box', 'on');
grid(axD, 'on'); grid(axD, 'minor'); box(axD, 'on');

axE = nexttile(tl1, [1 3]);
hold(axE, 'on');
bandIdx = simN > 0 & isfinite(simQ05) & isfinite(simQ95);
if any(bandIdx)
    bandX = [distBinCtr(bandIdx); flipud(distBinCtr(bandIdx))];
    bandY = [simQ05(bandIdx); flipud(simQ95(bandIdx))];
    patch(axE, 'XData', bandX, 'YData', bandY, 'FaceColor', [0.82 0.86 0.92], ...
        'FaceAlpha', 0.8, 'EdgeColor', 'none');
end
plot(axE, distBinCtr, simQ50, 'k-', 'LineWidth', 1.8);
compatTrue = Fd.IsEnvelopeCompatible;
compatFalse = ~Fd.IsEnvelopeCompatible;
scatter(axE, Fd.DistanceMid_m(compatFalse), Fd.Duration_ms(compatFalse), 10, [0.25 0.25 0.70], 'filled', ...
    'MarkerFaceAlpha', 0.55, 'MarkerEdgeAlpha', 0.55);
scatter(axE, Fd.DistanceMid_m(compatTrue), Fd.Duration_ms(compatTrue), 12, [0.82 0.72 0.10], 'filled', ...
    'MarkerFaceAlpha', 0.80, 'MarkerEdgeAlpha', 0.80);
hold(axE, 'off');
xlabel(axE, 'Inferred target-distance interval midpoint (m)', 'Interpreter', 'latex');
ylabel(axE, 'Call duration $T_c$ (ms)', 'Interpreter', 'latex');
title(axE, sprintf(['\\textbf{E. Envelope compatibility: ', ...
    '%.1f\\%% (IQR %.1f--%.1f\\%%)}'], ...
    compatPctMedian, compatPctQ1, compatPctQ3), ...
    'Interpreter', 'latex');
lgdEnvelope = legend(axE, {'sim. 5--95\% band', 'sim. median', 'incompatible', 'compatible'}, ...
    'Location', 'best', 'Box', 'off');
grid(axE, 'on'); grid(axE, 'minor'); box(axE, 'on');

formatLatex(fig1, "compact-landscape");
reserveLegendMetricMargin(lgdDuration, 0.10, axD);
reserveLegendMetricMargin(lgdEnvelope, 0.10, axE);

if savePlots
    exportPaperFigure(fig1, fullfile(outDir, 'field_vs_sim_observables'));
end
if saveStats
    writetable(velocityOverallSummary, fullfile(outDir, 'field_sim_effective_velocity_summary.csv'));
    writetable(velocityBandSummary, fullfile(outDir, 'field_sim_velocity_band_summary.csv'));
    writetable(vEffFieldSummary, fullfile(outDir, 'field_effective_velocity_binned_summary.csv'));
    writetable(vEffSimSummary, fullfile(outDir, 'simulation_effective_velocity_binned_summary.csv'));
    writetable(timingInteractionSummary, fullfile(outDir, 'field_sim_timing_interaction_summary.csv'));
    writetable(timingResidualSummary, fullfile(outDir, 'field_sim_timing_residual_summary.csv'));
    writetable(timingResidualBinSummary, fullfile(outDir, 'field_sim_timing_residual_binned_summary.csv'));
    writetable(durationCompatibilitySummary, fullfile(outDir, ...
        'duration_envelope_compatibility_seed_summary.csv'));
    writetable(table(durationCompatibilitySeeds(:), ...
        compatPctReplicates, 'VariableNames', ...
        {'SimulationSeed', 'CompatibleFieldCalls_pct'}), ...
        fullfile(outDir, ...
        'duration_envelope_compatibility_by_seed.csv'));
end
end

%% Lower-envelope test for duration versus rate
if makeDurationEnvelopeTestFigure
    figDurationEnvelope = figure('Color', 'w');
    tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    scatter(T.CallRate_Hz(T.IsKeep), T.Duration_ms(T.IsKeep), 10, T.Distance_m(T.IsKeep), 'filled', ...
        'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0.35);
    hold on;
    if showFlaggedPoints
        scatter(T.CallRate_Hz(T.IsAnyOutlier), T.Duration_ms(T.IsAnyOutlier), 14, ...
            [0.85 0.10 0.10], 'x', 'LineWidth', 0.8);
    end
    plotDurationFitCurve(gca, fieldDurationFit, fieldDurXLim, 'k-', 1.6);
    plotDurationFitCurve(gca, fieldDurationMedianFit, fieldDurXLim, 'b-.', 1.8);
    plotDurationFitCurve(gca, fieldDurationEnvelopeFit, fieldDurXLim, 'r--', 1.8);
    plotEnvelopeSupportPoints(fieldDurationMedianFit, 'b');
    plotEnvelopeSupportPoints(fieldDurationEnvelopeFit, 'r');
    hold off;
    xlabel('Call rate (Hz)', 'Interpreter', 'latex');
    ylabel('Call duration (ms)', 'Interpreter', 'latex');
    title('\textbf{A. Field: central, median, and lower-envelope fits}', 'Interpreter', 'latex');
    xlim(fieldDurXLim);
    ylim(fieldDurYLim);
    legend({'field data','central hyperbola','binwise-median hyperbola', ...
        'lower-envelope hyperbola','median bin points','envelope bin points'}, ...
        'Location', 'northeast');
    grid on; grid minor; box on;

    nexttile;
    scatter(S.CallRate_Hz, S.Tcall_ms, 10, S.InterCallDisplacement_m, 'filled', ...
        'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0.35);
    hold on;
    plotDurationFitCurve(gca, simDurationFit, fieldDurXLim, 'k-', 1.6);
    plotDurationFitCurve(gca, simDurationMedianFit, fieldDurXLim, 'b-.', 1.8);
    plotDurationFitCurve(gca, simDurationEnvelopeFit, fieldDurXLim, 'r--', 1.8);
    plotEnvelopeSupportPoints(simDurationMedianFit, 'b');
    plotEnvelopeSupportPoints(simDurationEnvelopeFit, 'r');
    hold off;
    xlabel('Call rate (Hz)', 'Interpreter', 'latex');
    ylabel('Call duration (ms)', 'Interpreter', 'latex');
    title('\textbf{B. Simulation: central, median, and lower-envelope fits}', 'Interpreter', 'latex');
    xlim(fieldDurXLim);
    ylim(fieldDurYLim);
    legend({'simulation data','central hyperbola','binwise-median hyperbola', ...
        'lower-envelope hyperbola','median bin points','envelope bin points'}, ...
        'Location', 'northeast');
    grid on; grid minor; box on;

    cb = colorbar;
    cb.Label.String = 'Distance or inter-call displacement';

    formatLatex(figDurationEnvelope, "full-landscape");

    if savePlots
        exportPaperFigure(figDurationEnvelope, fullfile(outDir, 'duration_rate_envelope_test'));
    end
end

%% Field versus simulation: IPI and timing decomposition
if makeLegacyTimingDecompositionFigure
fig2 = figure('Color', 'w', 'Position', [100 100 1000 600]);
tl2 = tiledlayout(fig2, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

fieldTaXLim = [min(T.Distance_m(T.IsKeep), [], 'omitnan') max(T.Distance_m(T.IsKeep), [], 'omitnan')];
fieldTaYLim = [min(T.Ta_fromIPI_ms(T.IsKeep), [], 'omitnan') max(T.Ta_fromIPI_ms(T.IsKeep), [], 'omitnan')];
fieldTbXLim = [min(T.Distance_m(T.IsKeep), [], 'omitnan') max(T.Distance_m(T.IsKeep), [], 'omitnan')];
fieldTbYLim = [min(T.Tb_fromIPI_ms(T.IsKeep), [], 'omitnan') max(T.Tb_fromIPI_ms(T.IsKeep), [], 'omitnan')];

ax2A = nexttile(tl2, 1);
scatter(ax2A, T.Distance_m(T.IsKeep), T.IPI_ms(T.IsKeep), 10, T.Velocity_m_s(T.IsKeep), 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75);
caxis(ax2A, velocityColourLimits_m_s);
if showFlaggedPoints
    hold(ax2A, 'on');
    scatter(ax2A, T.Distance_m(T.IsAnyOutlier), T.IPI_ms(T.IsAnyOutlier), 14, ...
        [0.85 0.10 0.10], 'x', 'LineWidth', 0.8);
end
% xlabel(ax2A, 'Inter-call displacement (m)', 'Interpreter', 'latex');
xlabel('');
ylabel(ax2A, 'IPI (ms)', 'Interpreter', 'latex');
title(ax2A, '\textbf{A. Field IPI}', 'Interpreter', 'latex');
hold(ax2A, 'on');
validIpi = T.IsKeep & isfinite(T.Distance_m) & isfinite(T.IPI_ms);
if sum(validIpi) >= 2
    pIpi = polyfit(T.Distance_m(validIpi), T.IPI_ms(validIpi), 1);
    xIpi = linspace(min(T.Distance_m(validIpi)), max(T.Distance_m(validIpi)), 200);
    yIpi = polyval(pIpi, xIpi);
    plot(ax2A, xIpi, yIpi, 'k-', 'LineWidth', 1.6);
end
hold(ax2A, 'off');
xlim(ax2A, fieldIpiXLim);
ylim(ax2A, fieldIpiYLim);
grid(ax2A, 'on'); grid(ax2A, 'minor'); box(ax2A, 'on');
addInlineVelocityScaleStrip(ax2A, velocityColourLimits_m_s, fieldIpiXLim, fieldIpiYLim, 'Velocity');

ax2B = nexttile(tl2, 2);
scatter(ax2B, T.Distance_m(T.IsKeep), T.Ta_fromIPI_ms(T.IsKeep), 10, T.Velocity_m_s(T.IsKeep), 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75);
caxis(ax2B, velocityColourLimits_m_s);
if showFlaggedPoints
    hold(ax2B, 'on');
    scatter(ax2B, T.Distance_m(T.IsAnyOutlier), T.Ta_fromIPI_ms(T.IsAnyOutlier), 14, ...
        [0.85 0.10 0.10], 'x', 'LineWidth', 0.8);
end
% xlabel(ax2B, 'Inter-call displacement (m)', 'Interpreter', 'latex');
xlabel('');
ylabel(ax2B, '$T_a$ from IPI (ms)', 'Interpreter', 'latex');
title(ax2B, sprintf('\\textbf{B. Field $T_a$ ($k_r = %g$)}', krModel), 'Interpreter', 'latex');
hold(ax2B, 'on');
validTa = T.IsKeep & isfinite(T.Distance_m) & isfinite(T.Ta_fromIPI_ms);
if sum(validTa) >= 2
    pTa = polyfit(T.Distance_m(validTa), T.Ta_fromIPI_ms(validTa), 1);
    xTa = linspace(min(T.Distance_m(validTa)), max(T.Distance_m(validTa)), 200);
    yTa = polyval(pTa, xTa);
    plot(ax2B, xTa, yTa, 'k-', 'LineWidth', 1.6);
end
hold(ax2B, 'off');
xlim(ax2B, fieldTaXLim);
ylim(ax2B, fieldTaYLim);
grid(ax2B, 'on'); grid(ax2B, 'minor'); box(ax2B, 'on');

ax2C = nexttile(tl2, 3);
scatter(ax2C, T.Distance_m(T.IsKeep), T.Tb_fromIPI_ms(T.IsKeep), 10, T.Velocity_m_s(T.IsKeep), 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75);
caxis(ax2C, velocityColourLimits_m_s);
hold(ax2C, 'on');
validTb = T.IsKeep & isfinite(T.Distance_m) & isfinite(T.Tb_fromIPI_ms);
if sum(validTb) >= 2
    pTb = polyfit(T.Distance_m(validTb), T.Tb_fromIPI_ms(validTb), 1);
    xTb = linspace(min(T.Distance_m(validTb)), max(T.Distance_m(validTb)), 200);
    yTb = polyval(pTb, xTb);
    plot(ax2C, xTb, yTb, 'k-', 'LineWidth', 1.6);
end
hold(ax2C, 'off');
% xlabel(ax2C, 'Inter-call displacement (m)', 'Interpreter', 'latex');
xlabel('');
ylabel(ax2C, '$T_b$ from IPI (ms)', 'Interpreter', 'latex');
title(ax2C, sprintf('\\textbf{C. Field $T_b$ ($k_r = %g$)}', krModel), 'Interpreter', 'latex');
xlim(ax2C, fieldTbXLim);
ylim(ax2C, fieldTbYLim);
grid(ax2C, 'on'); grid(ax2C, 'minor'); box(ax2C, 'on');

ax2D = nexttile(tl2, 4);
scatter(ax2D, S.InterCallDisplacement_m, S.IPI_ms, 10, S.BatSpeed_m_s, 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75);
caxis(ax2D, velocityColourLimits_m_s);
hold(ax2D, 'on');
validIpiSim = isfinite(S.InterCallDisplacement_m) & isfinite(S.IPI_ms);
if sum(validIpiSim) >= 2
    pIpiSim = polyfit(S.InterCallDisplacement_m(validIpiSim), S.IPI_ms(validIpiSim), 1);
    xIpiSim = linspace(min(S.InterCallDisplacement_m(validIpiSim)), max(S.InterCallDisplacement_m(validIpiSim)), 200);
    yIpiSim = polyval(pIpiSim, xIpiSim);
    plot(ax2D, xIpiSim, yIpiSim, 'k-', 'LineWidth', 1.6);
end
hold(ax2D, 'off');
xlabel(ax2D, 'Inter-call displacement (m)', 'Interpreter', 'latex');
ylabel(ax2D, 'IPI (ms)', 'Interpreter', 'latex');
title(ax2D, '\textbf{D. Simulated IPI}', 'Interpreter', 'latex');
xlim(ax2D, fieldIpiXLim);
ylim(ax2D, fieldIpiYLim);
grid(ax2D, 'on'); grid(ax2D, 'minor'); box(ax2D, 'on');
addInlineVelocityScaleStrip(ax2D, velocityColourLimits_m_s, fieldIpiXLim, fieldIpiYLim, 'Velocity');

ax2E = nexttile(tl2, 5);
scatter(ax2E, S.InterCallDisplacement_m, S.Ta_ms, 10, S.BatSpeed_m_s, 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75);
caxis(ax2E, velocityColourLimits_m_s);
hold(ax2E, 'on');
validTaSim = isfinite(S.InterCallDisplacement_m) & isfinite(S.Ta_ms);
if sum(validTaSim) >= 2
    pTaSim = polyfit(S.InterCallDisplacement_m(validTaSim), S.Ta_ms(validTaSim), 1);
    xTaSim = linspace(min(S.InterCallDisplacement_m(validTaSim)), max(S.InterCallDisplacement_m(validTaSim)), 200);
    yTaSim = polyval(pTaSim, xTaSim);
    plot(ax2E, xTaSim, yTaSim, 'k-', 'LineWidth', 1.6);
end
hold(ax2E, 'off');
xlabel(ax2E, 'Inter-call displacement (m)', 'Interpreter', 'latex');
ylabel(ax2E, '$T_a$ from simulation (ms)', 'Interpreter', 'latex');
title(ax2E, '\textbf{E. Simulated $T_a$}', 'Interpreter', 'latex');
xlim(ax2E, fieldTaXLim);
ylim(ax2E, fieldTaYLim);
grid(ax2E, 'on'); grid(ax2E, 'minor'); box(ax2E, 'on');

ax2F = nexttile(tl2, 6);
scatter(ax2F, S.InterCallDisplacement_m, S.Tb_ms, 10, S.BatSpeed_m_s, 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75);
caxis(ax2F, velocityColourLimits_m_s);
hold(ax2F, 'on');
validTbSim = isfinite(S.InterCallDisplacement_m) & isfinite(S.Tb_ms);
if sum(validTbSim) >= 2
    pTbSim = polyfit(S.InterCallDisplacement_m(validTbSim), S.Tb_ms(validTbSim), 1);
    xTbSim = linspace(min(S.InterCallDisplacement_m(validTbSim)), max(S.InterCallDisplacement_m(validTbSim)), 200);
    yTbSim = polyval(pTbSim, xTbSim);
    plot(ax2F, xTbSim, yTbSim, 'k-', 'LineWidth', 1.6);
end
hold(ax2F, 'off');
xlabel(ax2F, 'Inter-call displacement (m)', 'Interpreter', 'latex');
ylabel(ax2F, '$T_b$ from simulation (ms)', 'Interpreter', 'latex');
title(ax2F, '\textbf{F. Simulated $T_b$}', 'Interpreter', 'latex');
xlim(ax2F, fieldTbXLim);
ylim(ax2F, fieldTbYLim);
grid(ax2F, 'on'); grid(ax2F, 'minor'); box(ax2F, 'on');

formatLatex(fig2, "compact-landscape");

if savePlots
    exportPaperFigure(fig2, fullfile(outDir, 'field_vs_sim_timing_kr5'));
end
end

%% Independent observables and conditional framework inferences
% This exploratory figure separates directly observed/matched quantities
% from quantities inferred under a fixed k_r. In particular, field T_a and
% T_b are not treated as independent measurements: both are deterministic
% transforms of IPI when k_r is fixed.
runConditionalDiagnostics = makeTimingDiagnosticsFigure || ...
    makeIndependentDiagnosticFigure || makeExtendedPhiProfileFigure || ...
    makeD0PhiProfileFigure || makeDistanceBinnedEchoWindowFigure || ...
    runD0PhiSensitivity;
if runConditionalDiagnostics
    fieldDiagRows = T.IsKeep & isfinite(T.Distance_m) & T.Distance_m > 0 & ...
        isfinite(T.IPI_s) & T.IPI_s > 0 & ...
        isfinite(T.Velocity_m_s) & T.Velocity_m_s >= 0 & ...
        isfinite(T.Duration_s) & T.Duration_s > 0;
    simDiagRows = isfinite(S.InterCallDisplacement_m) & S.InterCallDisplacement_m > 0 & ...
        isfinite(S.IPI_s) & S.IPI_s > 0 & ...
        isfinite(S.BatSpeed_m_s) & S.BatSpeed_m_s >= 0 & ...
        isfinite(S.Tcall_s) & S.Tcall_s > 0 & ...
        isfinite(S.AnchorDistance_m) & S.AnchorDistance_m >= 0 & ...
        isfinite(S.RelativeClosingVelocityForTiming_m_s);

    fieldDiag = T(fieldDiagRows, :);
    simDiag = S(simDiagRows, :);

    % Shared IPI-displacement summaries.
    diagDispMax = max([fieldDiag.Distance_m; simDiag.InterCallDisplacement_m], [], 'omitnan');
    diagDispEdges = linspace(0, diagDispMax, 18);
    [fieldDispCtr, fieldIpiQ25, fieldIpiQ50, fieldIpiQ75] = ...
        binnedQuantiles(fieldDiag.Distance_m, fieldDiag.IPI_ms, diagDispEdges, [25 50 75]);
    [simDispCtr, simIpiQ25, simIpiQ50, simIpiQ75] = ...
        binnedQuantiles(simDiag.InterCallDisplacement_m, simDiag.IPI_ms, diagDispEdges, [25 50 75]);

    % Shared velocity-conditioned IPI summaries. These remain independent
    % of the fixed-k_r timing decomposition.
    diagVelocityEdges = velocityColourLimits_m_s(1):velocityMatchBinWidth_m_s: ...
        velocityColourLimits_m_s(2);
    [fieldVelocityCtr, fieldVelocityIpiQ25, fieldVelocityIpiQ50, fieldVelocityIpiQ75] = ...
        binnedQuantiles(fieldDiag.Velocity_m_s, fieldDiag.IPI_ms, ...
        diagVelocityEdges, [25 50 75]);
    [simVelocityCtr, simVelocityIpiQ25, simVelocityIpiQ50, simVelocityIpiQ75] = ...
        binnedQuantiles(simDiag.BatSpeed_m_s, simDiag.IPI_ms, ...
        diagVelocityEdges, [25 50 75]);
    fieldVelocityBinCount = histcounts(fieldDiag.Velocity_m_s, diagVelocityEdges)';
    simVelocityBinCount = histcounts(simDiag.BatSpeed_m_s, diagVelocityEdges)';

    % Infer the acquisition window from the observable IPI under k_r.
    simTaInferred_s = simDiag.IPI_s ./ (1 + krModel);
    simRangeFactor = 0.5 .* (cSound + simDiag.RelativeClosingVelocityForTiming_m_s);
    simDistanceTrue_m = simDiag.AnchorDistance_m;

    % With phi unknown, phi in [0,1] defines a near--far distance interval.
    simDistanceFar_m = simRangeFactor .* simTaInferred_s;
    simDistanceNear_m = simRangeFactor .* max(simTaInferred_s - simDiag.Tcall_s, 0);
    coverageTolerance_m = 1e-9;
    simDistanceCoverage = simDistanceTrue_m >= simDistanceNear_m - coverageTolerance_m & ...
        simDistanceTrue_m <= simDistanceFar_m + coverageTolerance_m;
    simDistanceCoveragePct = 100 * mean(simDistanceCoverage, 'omitnan');
    simDistanceIntervalWidth_m = simDistanceFar_m - simDistanceNear_m;
    simDistanceMedianIntervalWidth_m = median(simDistanceIntervalWidth_m, 'omitnan');
    distanceUncertaintyEdges = linspace(0, prctile(simDistanceTrue_m, 99.5), 14);
    [distanceUncertaintyCtr, distanceWidthQ25, distanceWidthQ50, distanceWidthQ75] = ...
        binnedQuantiles(simDistanceTrue_m, simDistanceIntervalWidth_m, ...
        distanceUncertaintyEdges, [25 50 75]);

    % Profile phi as a generative parameter. At every candidate phi, build
    % a fresh simulation ensemble, match it to the field velocity
    % distribution, and compare distributions using the mean absolute gap
    % between their 1st--99th percentiles.
    phiProfileSimulationConfig = struct( ...
        'CacheVersion', 1, ...
        'Kr', krModel, ...
        'MaxVelocity_m_s', maxVelocityForKeep_m_s, ...
        'PhiValues', phiProfileValues, ...
        'NumSequencesPerPhi', phiProfileNumSequencesPerPhi, ...
        'InitialSpeedRange_m_s', simulationInitialSpeedRange_m_s, ...
        'SeedBase', phiProfileSimulationSeedBase, ...
        'InitialDistance_m', defaultSimulationInitialDistance_m);
    phiProfileSimulationCacheFile = fullfile(simulationCacheDir, ...
        'phi_profile_main.mat');
    [SphiProfileAll, phiProfileSimulationLoaded] = ...
        loadSimulationCache(phiProfileSimulationCacheFile, ...
        phiProfileSimulationConfig, loadSavedSimulationData);
    if ~phiProfileSimulationLoaded
        SphiProfileAll = runFieldValidationSimulation(thisDir, krModel, ...
            maxVelocityForKeep_m_s, phiProfileValues, ...
            phiProfileNumSequencesPerPhi, ...
            simulationInitialSpeedRange_m_s, ...
            phiProfileSimulationSeedBase, ...
            defaultSimulationInitialDistance_m);
    end
    phiProfileNeedsCacheWrite = ~phiProfileSimulationLoaded || ...
        width(SphiProfileAll) > numel(fieldValidationCacheVariables());
    SphiProfileAll = ...
        trimFieldValidationSimulationTable(SphiProfileAll);
    if phiProfileNeedsCacheWrite
        saveSimulationCache(phiProfileSimulationCacheFile, ...
            SphiProfileAll, phiProfileSimulationConfig, ...
            saveSimulationData);
    end
    fprintf('Phi-profile simulation source: %s\n', ...
        simulationSourceLabel(phiProfileSimulationLoaded));
    numPhiProfile = numel(phiProfileValues);
    profileVelocityMatchEdges_m_s = profileVelocityLimits_m_s(1): ...
        profileVelocityMatchBinWidth_m_s:profileVelocityLimits_m_s(2);
    numProfileVelocityBins = numel(profileVelocityMatchEdges_m_s) - 1;
    fieldProfileBinCount = histcounts(fieldDiag.Velocity_m_s, ...
        profileVelocityMatchEdges_m_s)';
    fieldProfileBinCount(fieldProfileBinCount < ...
        diagnosticMinVelocityBinCount) = 0;
    phiAvailableBinCount = zeros(numProfileVelocityBins, numPhiProfile);
    for p = 1:numPhiProfile
        phiRows = abs(SphiProfileAll.Phi - phiProfileValues(p)) < 1e-12;
        phiAvailableBinCount(:, p) = histcounts( ...
            SphiProfileAll.BatSpeed_m_s(phiRows), ...
            profileVelocityMatchEdges_m_s)';
    end
    [profileTargetBinCount, phiProfileFieldFraction] = ...
        scaledVelocityTargetCounts(fieldProfileBinCount, ...
        phiAvailableBinCount);
    [fieldProfileRows, ~] = ...
        sampleRowsByVelocityTarget(fieldDiag.Velocity_m_s, ...
        profileVelocityMatchEdges_m_s, profileTargetBinCount, velocityMatchSeed);
    fieldProfile = fieldDiag(fieldProfileRows, :);
    fieldProfileTa_s = fieldProfile.IPI_s ./ (1 + krModel);
    fieldProfileRangeFactor = 0.5 .* (cSound + fieldProfile.Velocity_m_s);
    fieldProfileTb_ms = 1000 .* (krModel ./ (1 + krModel)) .* fieldProfile.IPI_s;

    distancePhiLoss_m = nan(numPhiProfile, 1);
    tbPhiLoss_ms = nan(numPhiProfile, 1);
    phiProfileMatchedCount = zeros(numPhiProfile, 1);
    phiProfileShortage = zeros(numPhiProfile, 1);
    fieldDistanceMedian_m = nan(numPhiProfile, 1);
    simDistanceMedian_m = nan(numPhiProfile, 1);
    simTbMedian_ms = nan(numPhiProfile, 1);
    fieldDistanceClippedPct = nan(numPhiProfile, 1);

    for p = 1:numPhiProfile
        phiCandidate = phiProfileValues(p);
        phiRows = abs(SphiProfileAll.Phi - phiCandidate) < 1e-12;
        SphiRaw = SphiProfileAll(phiRows, :);
        [simProfileRows, simProfileActualBinCount] = ...
            sampleRowsByVelocityTarget(SphiRaw.BatSpeed_m_s, ...
            profileVelocityMatchEdges_m_s, profileTargetBinCount, ...
            velocityMatchSeed);
        SphiMatched = SphiRaw(simProfileRows, :);

        phiProfileMatchedCount(p) = height(SphiMatched);
        phiProfileShortage(p) = sum( ...
            profileTargetBinCount - simProfileActualBinCount);

        fieldDelayCandidate_s = fieldProfileTa_s - ...
            phiCandidate .* fieldProfile.Duration_s;
        fieldDistanceClippedPct(p) = 100 .* mean(fieldDelayCandidate_s <= 0);
        fieldDistanceCandidate_m = fieldProfileRangeFactor .* ...
            max(fieldDelayCandidate_s, 0);
        simDistanceCandidate_m = SphiMatched.AnchorDistance_m;
        simTbCandidate_ms = 1000 .* SphiMatched.TbEffective_s;

        distancePhiLoss_m(p) = meanAbsoluteQuantileGap( ...
            fieldDistanceCandidate_m, simDistanceCandidate_m, phiProfileQuantiles);
        tbPhiLoss_ms(p) = meanAbsoluteQuantileGap( ...
            fieldProfileTb_ms, simTbCandidate_ms, phiProfileQuantiles);
        fieldDistanceMedian_m(p) = median(fieldDistanceCandidate_m, 'omitnan');
        simDistanceMedian_m(p) = median(simDistanceCandidate_m, 'omitnan');
        simTbMedian_ms(p) = median(simTbCandidate_ms, 'omitnan');
    end

    [bestDistanceLoss_m, bestDistanceIdx] = min(distancePhiLoss_m);
    [bestTbLoss_ms, bestTbIdx] = min(tbPhiLoss_ms);
    bestDistancePhi = phiProfileValues(bestDistanceIdx);
    bestTbPhi = phiProfileValues(bestTbIdx);

    % At the distance-profile optimum, pair field and simulated calls by
    % acquisition-time and anchor-distance rank within each velocity bin.
    % The resulting required phi is a conditional distribution-matching
    % diagnostic, not a directly observed call-specific estimate.
    bestPhiRows = abs(SphiProfileAll.Phi - bestDistancePhi) < 1e-12;
    SbestPhiRaw = SphiProfileAll(bestPhiRows, :);
    [bestPhiSampleRows, ~] = sampleRowsByVelocityTarget( ...
        SbestPhiRaw.BatSpeed_m_s, profileVelocityMatchEdges_m_s, ...
        profileTargetBinCount, velocityMatchSeed);
    SbestPhi = SbestPhiRaw(bestPhiSampleRows, :);
    [requiredPhi, pairedAnchorDistance_m, pairedCallDuration_ms, ...
        pairedTa_ms, pairedDelay_ms] = ...
        requiredPhiByVelocityRankMatching(fieldProfile, SbestPhi, ...
        profileVelocityMatchEdges_m_s, krModel, cSound);
    pairedPhiWindow_ms = requiredPhi .* pairedCallDuration_ms;

    phiDistanceMax_m = ceil(prctile(pairedAnchorDistance_m, 99) ./ ...
        phiDistanceBinWidth_m) .* phiDistanceBinWidth_m;
    phiDistanceEdges_m = 0:phiDistanceBinWidth_m:phiDistanceMax_m;
    [phiDistanceCtr_m, requiredPhiQ25, requiredPhiQ50, requiredPhiQ75] = ...
        binnedQuantiles(pairedAnchorDistance_m, requiredPhi, ...
        phiDistanceEdges_m, [25 50 75]);
    [~, callDurationQ25_ms, callDurationQ50_ms, callDurationQ75_ms] = ...
        binnedQuantiles(pairedAnchorDistance_m, pairedCallDuration_ms, ...
        phiDistanceEdges_m, [25 50 75]);
    [~, phiWindowQ25_ms, phiWindowQ50_ms, phiWindowQ75_ms] = ...
        binnedQuantiles(pairedAnchorDistance_m, pairedPhiWindow_ms, ...
        phiDistanceEdges_m, [25 50 75]);
    [~, pairedDelayQ25_ms, pairedDelayQ50_ms, pairedDelayQ75_ms] = ...
        binnedQuantiles(pairedAnchorDistance_m, pairedDelay_ms, ...
        phiDistanceEdges_m, [25 50 75]);
    [~, pairedTaQ25_ms, pairedTaQ50_ms, pairedTaQ75_ms] = ...
        binnedQuantiles(pairedAnchorDistance_m, pairedTa_ms, ...
        phiDistanceEdges_m, [25 50 75]);
    requiredPhiBinCount = histcounts(pairedAnchorDistance_m, ...
        phiDistanceEdges_m)';
    requiredPhiBinKeep = isfinite(phiDistanceCtr_m) & ...
        isfinite(requiredPhiQ25) & isfinite(requiredPhiQ50) & ...
        isfinite(requiredPhiQ75) & ...
        requiredPhiBinCount >= diagnosticMinVelocityBinCount;
    absoluteTimingBinKeep = requiredPhiBinKeep & ...
        isfinite(callDurationQ25_ms) & isfinite(callDurationQ50_ms) & ...
        isfinite(callDurationQ75_ms) & ...
        isfinite(phiWindowQ25_ms) & isfinite(phiWindowQ50_ms) & ...
        isfinite(phiWindowQ75_ms) & ...
        isfinite(pairedDelayQ50_ms) & isfinite(pairedTaQ50_ms);

    nearTimingBinIdx = find(absoluteTimingBinKeep, 1, 'first');
    farTimingBinIdx = find(absoluteTimingBinKeep, 1, 'last');

    % Two-parameter sensitivity: profile phi independently at each initial
    % simulation distance using common random numbers and one common
    % velocity-stratified support across the complete d0-by-phi grid.
    if runD0PhiSensitivity
        numD0Profile = numel(d0ProfileValues_m);
        numD0PhiProfile = numel(d0PhiProfileValues);
        d0PhiSimulationConfig = struct( ...
            'CacheVersion', 1, ...
            'Kr', krModel, ...
            'MaxVelocity_m_s', maxVelocityForKeep_m_s, ...
            'D0Values_m', d0ProfileValues_m, ...
            'PhiValues', d0PhiProfileValues, ...
            'NumSequencesPerCombination', ...
                d0PhiNumSequencesPerCombination, ...
            'InitialSpeedRange_m_s', simulationInitialSpeedRange_m_s, ...
            'SeedBase', d0PhiSimulationSeedBase);
        d0PhiSimulationCacheFile = fullfile(simulationCacheDir, ...
            'd0_phi_profile.mat');
        [d0PhiSimulationTables, d0PhiSimulationLoaded] = ...
            loadSimulationCache(d0PhiSimulationCacheFile, ...
            d0PhiSimulationConfig, loadSavedSimulationData);
        if ~d0PhiSimulationLoaded
            d0PhiSimulationTables = cell(numD0Profile, 1);
        end
        d0PhiSimulationNeedsCacheWrite = ~d0PhiSimulationLoaded;
        d0PhiAvailableBinCount = zeros(numProfileVelocityBins, ...
            numD0Profile, numD0PhiProfile);

        fprintf('\nd0-by-phi simulation source: %s\n', ...
            simulationSourceLabel(d0PhiSimulationLoaded));
        if ~d0PhiSimulationLoaded
            fprintf(['Running d0-by-phi profile: %d starting distances ', ...
                'x %d phi values\n'], numD0Profile, numD0PhiProfile);
        end
        for d = 1:numD0Profile
            if ~d0PhiSimulationLoaded
                d0PhiSimulationTables{d} = ...
                    runFieldValidationSimulation( ...
                    thisDir, krModel, maxVelocityForKeep_m_s, ...
                    d0PhiProfileValues, ...
                    d0PhiNumSequencesPerCombination, ...
                    simulationInitialSpeedRange_m_s, ...
                    d0PhiSimulationSeedBase, d0ProfileValues_m(d));
            end
            d0PhiSimulationNeedsCacheWrite = ...
                d0PhiSimulationNeedsCacheWrite || ...
                width(d0PhiSimulationTables{d}) > ...
                numel(fieldValidationCacheVariables());
            d0PhiSimulationTables{d} = ...
                trimFieldValidationSimulationTable( ...
                d0PhiSimulationTables{d});
            for p = 1:numD0PhiProfile
                phiRows = abs(d0PhiSimulationTables{d}.Phi - ...
                    d0PhiProfileValues(p)) < 1e-12;
                d0PhiAvailableBinCount(:, d, p) = histcounts( ...
                    d0PhiSimulationTables{d}.BatSpeed_m_s(phiRows), ...
                    profileVelocityMatchEdges_m_s)';
            end
        end
        if d0PhiSimulationNeedsCacheWrite
            saveSimulationCache(d0PhiSimulationCacheFile, ...
                d0PhiSimulationTables, d0PhiSimulationConfig, ...
                saveSimulationData);
        end

        d0PhiAvailabilityFlat = reshape(d0PhiAvailableBinCount, ...
            numProfileVelocityBins, []);
        [d0PhiTargetBinCount, d0PhiFieldFraction] = ...
            scaledVelocityTargetCounts(fieldProfileBinCount, ...
            d0PhiAvailabilityFlat);
        [d0PhiFieldRows, ~] = sampleRowsByVelocityTarget( ...
            fieldDiag.Velocity_m_s, profileVelocityMatchEdges_m_s, ...
            d0PhiTargetBinCount, velocityMatchSeed);
        d0PhiField = fieldDiag(d0PhiFieldRows, :);
        d0PhiFieldTa_s = d0PhiField.IPI_s ./ (1 + krModel);
        d0PhiFieldRangeFactor = 0.5 .* ...
            (cSound + d0PhiField.Velocity_m_s);
        d0PhiFieldTb_ms = 1000 .* (krModel ./ (1 + krModel)) .* ...
            d0PhiField.IPI_s;

        d0DistanceLoss_m = nan(numD0Profile, numD0PhiProfile);
        d0TbLoss_ms = nan(numD0Profile, numD0PhiProfile);
        d0PhiMatchedCount = zeros(numD0Profile, numD0PhiProfile);

        for d = 1:numD0Profile
            for p = 1:numD0PhiProfile
                phiCandidate = d0PhiProfileValues(p);
                phiRows = abs(d0PhiSimulationTables{d}.Phi - ...
                    phiCandidate) < 1e-12;
                SgridRaw = d0PhiSimulationTables{d}(phiRows, :);
                [gridSampleRows, ~] = sampleRowsByVelocityTarget( ...
                    SgridRaw.BatSpeed_m_s, profileVelocityMatchEdges_m_s, ...
                    d0PhiTargetBinCount, velocityMatchSeed);
                SgridMatched = SgridRaw(gridSampleRows, :);
                d0PhiMatchedCount(d, p) = height(SgridMatched);

                fieldDelayCandidate_s = d0PhiFieldTa_s - ...
                    phiCandidate .* d0PhiField.Duration_s;
                fieldDistanceCandidate_m = d0PhiFieldRangeFactor .* ...
                    max(fieldDelayCandidate_s, 0);
                d0DistanceLoss_m(d, p) = meanAbsoluteQuantileGap( ...
                    fieldDistanceCandidate_m, SgridMatched.AnchorDistance_m, ...
                    phiProfileQuantiles);
                d0TbLoss_ms(d, p) = meanAbsoluteQuantileGap( ...
                    d0PhiFieldTb_ms, 1000 .* SgridMatched.TbEffective_s, ...
                    phiProfileQuantiles);
            end
        end

        [bestDistanceLossByD0_m, bestDistancePhiIdxByD0] = ...
            min(d0DistanceLoss_m, [], 2);
        [bestTbLossByD0_ms, bestTbPhiIdxByD0] = ...
            min(d0TbLoss_ms, [], 2);
        bestDistancePhiByD0 = d0PhiProfileValues(bestDistancePhiIdxByD0)';
        bestTbPhiByD0 = d0PhiProfileValues(bestTbPhiIdxByD0)';

        [globalD0DistanceLoss_m, globalDistanceLinearIdx] = ...
            min(d0DistanceLoss_m(:));
        [globalDistanceD0Idx, globalDistancePhiIdx] = ...
            ind2sub(size(d0DistanceLoss_m), globalDistanceLinearIdx);
        [globalD0TbLoss_ms, globalTbLinearIdx] = ...
            min(d0TbLoss_ms(:));
        [globalTbD0Idx, globalTbPhiIdx] = ...
            ind2sub(size(d0TbLoss_ms), globalTbLinearIdx);
    end

    fprintf('\nIndependent-observable and conditional-inference diagnostics\n');
    fprintf('Simulation true-distance coverage by unknown-phi interval: %.2f%%\n', ...
        simDistanceCoveragePct);
    fprintf('Simulation median unknown-phi interval width: %.3f m\n', ...
        simDistanceMedianIntervalWidth_m);
    fprintf('Distance-profile best phi: %.2f (mean quantile gap %.3f m)\n', ...
        bestDistancePhi, bestDistanceLoss_m);
    fprintf('At distance-profile optimum, median distance: field %.3f m; simulation %.3f m\n', ...
        fieldDistanceMedian_m(bestDistanceIdx), simDistanceMedian_m(bestDistanceIdx));
    fprintf('Field zero-distance clipping at distance-profile optimum: %.2f%%\n', ...
        fieldDistanceClippedPct(bestDistanceIdx));
    fprintf('T_b-profile best generative phi: %.2f (mean quantile gap %.2f ms)\n', ...
        bestTbPhi, bestTbLoss_ms);
    fprintf('At T_b-profile optimum, median T_b: field %.2f ms; simulation %.2f ms\n', ...
        median(fieldProfileTb_ms, 'omitnan'), simTbMedian_ms(bestTbIdx));
    fprintf(['Phi-profile common velocity-matched calls per value: %d--%d; ', ...
        'field calls: %d; retained field fraction: %.3f; total shortages: %d\n'], ...
        min(phiProfileMatchedCount), max(phiProfileMatchedCount), ...
        height(fieldProfile), phiProfileFieldFraction, ...
        sum(phiProfileShortage));
    fprintf('Rank-matched required phi: median %.2f, IQR %.2f--%.2f\n', ...
        median(requiredPhi, 'omitnan'), prctile(requiredPhi, 25), ...
        prctile(requiredPhi, 75));
    if ~isempty(nearTimingBinIdx) && ~isempty(farTimingBinIdx)
        fprintf(['Far-to-near absolute timing: distance %.2f to %.2f m; ', ...
            'median phi %.2f to %.2f; T_c %.2f to %.2f ms; ', ...
            'phi*T_c %.2f to %.2f ms; T_delay %.2f to %.2f ms; ', ...
            'T_a %.2f to %.2f ms\n'], ...
            phiDistanceCtr_m(farTimingBinIdx), ...
            phiDistanceCtr_m(nearTimingBinIdx), ...
            requiredPhiQ50(farTimingBinIdx), ...
            requiredPhiQ50(nearTimingBinIdx), ...
            callDurationQ50_ms(farTimingBinIdx), ...
            callDurationQ50_ms(nearTimingBinIdx), ...
            phiWindowQ50_ms(farTimingBinIdx), ...
            phiWindowQ50_ms(nearTimingBinIdx), ...
            pairedDelayQ50_ms(farTimingBinIdx), ...
            pairedDelayQ50_ms(nearTimingBinIdx), ...
            pairedTaQ50_ms(farTimingBinIdx), ...
            pairedTaQ50_ms(nearTimingBinIdx));
    end
    if runD0PhiSensitivity
        fprintf('d0-by-phi common matched calls per combination: %d\n', ...
            min(d0PhiMatchedCount(:)));
        fprintf('d0-by-phi retained field fraction: %.3f\n', ...
            d0PhiFieldFraction);
        fprintf('Global distance optimum: d0 = %.1f m, phi = %.2f, gap = %.3f m\n', ...
            d0ProfileValues_m(globalDistanceD0Idx), ...
            d0PhiProfileValues(globalDistancePhiIdx), ...
            globalD0DistanceLoss_m);
        fprintf('Global T_b optimum: d0 = %.1f m, phi = %.2f, gap = %.2f ms\n', ...
            d0ProfileValues_m(globalTbD0Idx), ...
            d0PhiProfileValues(globalTbPhiIdx), globalD0TbLoss_ms);
        disp(table(d0ProfileValues_m(:), bestDistancePhiByD0, ...
            bestDistanceLossByD0_m, bestTbPhiByD0, bestTbLossByD0_ms, ...
            'VariableNames', {'d0_m', 'BestDistancePhi', ...
            'DistanceGap_m', 'BestTbPhi', 'TbGap_ms'}));
    end

    intervalColour = [0.30 0.62 0.52];

    if makeIndependentDiagnosticFigure
    fig3 = figure('Color', 'w', 'Position', [100 100 1000 600]);
    tl3 = tiledlayout(fig3, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    % A. Observable IPI structure as binned median and IQR.
    ax3A = nexttile(tl3, 1);
    hold(ax3A, 'on');
    fieldBandKeep = isfinite(fieldDispCtr) & isfinite(fieldIpiQ25) & ...
        isfinite(fieldIpiQ50) & isfinite(fieldIpiQ75);
    simBandKeep = isfinite(simDispCtr) & isfinite(simIpiQ25) & ...
        isfinite(simIpiQ50) & isfinite(simIpiQ75);
    fill(ax3A, [fieldDispCtr(fieldBandKeep); flipud(fieldDispCtr(fieldBandKeep))], ...
        [fieldIpiQ25(fieldBandKeep); flipud(fieldIpiQ75(fieldBandKeep))], ...
        fieldColour, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    fill(ax3A, [simDispCtr(simBandKeep); flipud(simDispCtr(simBandKeep))], ...
        [simIpiQ25(simBandKeep); flipud(simIpiQ75(simBandKeep))], ...
        simColour, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    hFieldIpi = plot(ax3A, fieldDispCtr(fieldBandKeep), fieldIpiQ50(fieldBandKeep), ...
        '-', 'Color', fieldColour, 'LineWidth', 1.8);
    hSimIpi = plot(ax3A, simDispCtr(simBandKeep), simIpiQ50(simBandKeep), ...
        '-', 'Color', simColour, 'LineWidth', 1.8);
    hold(ax3A, 'off');
    xlabel(ax3A, 'Inter-call displacement (m)', 'Interpreter', 'latex');
    ylabel(ax3A, 'IPI (ms)', 'Interpreter', 'latex');
    title(ax3A, '\textbf{A. Observable IPI structure}', 'Interpreter', 'latex');
    grid(ax3A, 'on'); grid(ax3A, 'minor'); box(ax3A, 'on');

    % B. Kinematic support shown independently of the k_r decomposition.
    ax3B = nexttile(tl3, 2);
    fieldVelocitySorted = sort(fieldDiag.Velocity_m_s);
    simVelocitySorted = sort(simDiag.BatSpeed_m_s);
    fieldVelocityCdf = (1:numel(fieldVelocitySorted))' ./ numel(fieldVelocitySorted);
    simVelocityCdf = (1:numel(simVelocitySorted))' ./ numel(simVelocitySorted);
    stairs(ax3B, fieldVelocitySorted, fieldVelocityCdf, ...
        'Color', fieldColour, 'LineWidth', 1.8);
    hold(ax3B, 'on');
    stairs(ax3B, simVelocitySorted, simVelocityCdf, ...
        'Color', simColour, 'LineWidth', 1.8);
    hold(ax3B, 'off');
    xlabel(ax3B, 'Velocity (m s$^{-1}$)', 'Interpreter', 'latex');
    ylabel(ax3B, 'Empirical cumulative probability', 'Interpreter', 'latex');
    title(ax3B, '\textbf{B. Velocity matching}', 'Interpreter', 'latex');
    xlim(ax3B, velocityColourLimits_m_s);
    ylim(ax3B, [0 1]);
    grid(ax3B, 'on'); grid(ax3B, 'minor'); box(ax3B, 'on');

    % C. Observable IPI structure conditioned on velocity.
    ax3C = nexttile(tl3, 3);
    hold(ax3C, 'on');
    fieldVelocityBandKeep = isfinite(fieldVelocityCtr) & ...
        isfinite(fieldVelocityIpiQ25) & isfinite(fieldVelocityIpiQ50) & ...
        isfinite(fieldVelocityIpiQ75) & ...
        fieldVelocityBinCount >= diagnosticMinVelocityBinCount;
    simVelocityBandKeep = isfinite(simVelocityCtr) & ...
        isfinite(simVelocityIpiQ25) & isfinite(simVelocityIpiQ50) & ...
        isfinite(simVelocityIpiQ75) & ...
        simVelocityBinCount >= diagnosticMinVelocityBinCount;
    fill(ax3C, ...
        [fieldVelocityCtr(fieldVelocityBandKeep); flipud(fieldVelocityCtr(fieldVelocityBandKeep))], ...
        [fieldVelocityIpiQ25(fieldVelocityBandKeep); flipud(fieldVelocityIpiQ75(fieldVelocityBandKeep))], ...
        fieldColour, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    fill(ax3C, ...
        [simVelocityCtr(simVelocityBandKeep); flipud(simVelocityCtr(simVelocityBandKeep))], ...
        [simVelocityIpiQ25(simVelocityBandKeep); flipud(simVelocityIpiQ75(simVelocityBandKeep))], ...
        simColour, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    plot(ax3C, fieldVelocityCtr(fieldVelocityBandKeep), ...
        fieldVelocityIpiQ50(fieldVelocityBandKeep), '-', ...
        'Color', fieldColour, 'LineWidth', 1.8);
    plot(ax3C, simVelocityCtr(simVelocityBandKeep), ...
        simVelocityIpiQ50(simVelocityBandKeep), '-', ...
        'Color', simColour, 'LineWidth', 1.8);
    hold(ax3C, 'off');
    xlabel(ax3C, 'Velocity (m s$^{-1}$)', 'Interpreter', 'latex');
    ylabel(ax3C, 'IPI (ms)', 'Interpreter', 'latex');
    title(ax3C, '\textbf{C. Velocity-conditioned IPI}', 'Interpreter', 'latex');
    xlim(ax3C, velocityColourLimits_m_s);
    grid(ax3C, 'on'); grid(ax3C, 'minor'); box(ax3C, 'on');

    % D. Distance uncertainty introduced when phi is not known.
    ax3D = nexttile(tl3, 4);
    uncertaintyKeep = isfinite(distanceUncertaintyCtr) & isfinite(distanceWidthQ25) & ...
        isfinite(distanceWidthQ50) & isfinite(distanceWidthQ75);
    fill(ax3D, ...
        [distanceUncertaintyCtr(uncertaintyKeep); flipud(distanceUncertaintyCtr(uncertaintyKeep))], ...
        [distanceWidthQ25(uncertaintyKeep); flipud(distanceWidthQ75(uncertaintyKeep))], ...
        intervalColour, 'FaceAlpha', 0.24, 'EdgeColor', 'none');
    hold(ax3D, 'on');
    plot(ax3D, distanceUncertaintyCtr(uncertaintyKeep), ...
        distanceWidthQ50(uncertaintyKeep), '-', ...
        'Color', intervalColour, 'LineWidth', 1.8);
    hold(ax3D, 'off');
    xlabel(ax3D, 'True simulated anchor distance (m)', 'Interpreter', 'latex');
    ylabel(ax3D, 'Near--far interval width (m)', 'Interpreter', 'latex');
    title(ax3D, '\textbf{D. Unknown-$\phi$ uncertainty}', 'Interpreter', 'latex');
    text(ax3D, 0.04, 0.94, ...
        sprintf('coverage $=%.1f$ percent', simDistanceCoveragePct), ...
        'Units', 'normalized', 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    text(ax3D, 0.04, 0.84, ...
        sprintf('median width $=%.3f$ m', simDistanceMedianIntervalWidth_m), ...
        'Units', 'normalized', 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    grid(ax3D, 'on'); grid(ax3D, 'minor'); box(ax3D, 'on');

    % E. Phi profile for distance-distribution agreement.
    ax3E = nexttile(tl3, 5);
    plot(ax3E, phiProfileValues, distancePhiLoss_m, '-', ...
        'Color', intervalColour, 'LineWidth', 1.8);
    hold(ax3E, 'on');
    plot(ax3E, bestDistancePhi, bestDistanceLoss_m, 'o', ...
        'Color', intervalColour, 'MarkerFaceColor', intervalColour, ...
        'MarkerSize', 6, 'LineWidth', 1.0);
    hold(ax3E, 'off');
    xlabel(ax3E, '$\phi$', 'Interpreter', 'latex');
    ylabel(ax3E, 'Mean quantile gap (m)', 'Interpreter', 'latex');
    title(ax3E, '\textbf{E. Distance fit across $\phi$}', ...
        'Interpreter', 'latex');
    xlim(ax3E, [phiProfileValues(1) phiProfileValues(end)]);
    addPhiProfileRangeCue(ax3E, phiOriginalSimulationLimit, ...
        phiProfileValues(end));
    text(ax3E, 0.04, 0.94, ...
        sprintf('minimum at $\\phi=%.2f$', bestDistancePhi), ...
        'Units', 'normalized', 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    text(ax3E, 0.04, 0.74, ...
        sprintf('zero-clipped $=%.1f$ percent', ...
        fieldDistanceClippedPct(bestDistanceIdx)), ...
        'Units', 'normalized', 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    text(ax3E, 0.04, 0.84, ...
        sprintf('minimum gap $=%.3f$ m', bestDistanceLoss_m), ...
        'Units', 'normalized', 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    grid(ax3E, 'on'); grid(ax3E, 'minor'); box(ax3E, 'on');

    % F. Phi profile for realised simulation T_b against field-inferred T_b.
    ax3F = nexttile(tl3, 6);
    plot(ax3F, phiProfileValues, tbPhiLoss_ms, '-', ...
        'Color', simColour, 'LineWidth', 1.8);
    hold(ax3F, 'on');
    plot(ax3F, bestTbPhi, bestTbLoss_ms, 'o', ...
        'Color', simColour, 'MarkerFaceColor', simColour, ...
        'MarkerSize', 6, 'LineWidth', 1.0);
    hold(ax3F, 'off');
    xlabel(ax3F, 'Generative $\phi$', 'Interpreter', 'latex');
    ylabel(ax3F, 'Mean quantile gap (ms)', 'Interpreter', 'latex');
    title(ax3F, '\textbf{F. $T_b$ fit across $\phi$}', ...
        'Interpreter', 'latex');
    xlim(ax3F, [phiProfileValues(1) phiProfileValues(end)]);
    addPhiProfileRangeCue(ax3F, phiOriginalSimulationLimit, ...
        phiProfileValues(end));
    text(ax3F, 0.04, 0.94, ...
        sprintf('minimum at $\\phi=%.2f$', bestTbPhi), ...
        'Units', 'normalized', 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    text(ax3F, 0.04, 0.84, ...
        sprintf('minimum gap $=%.2f$ ms', bestTbLoss_ms), ...
        'Units', 'normalized', 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    grid(ax3F, 'on'); grid(ax3F, 'minor'); box(ax3F, 'on');

    lgd3 = legend(ax3A, [hFieldIpi hSimIpi], {'field', 'simulation'}, ...
        'Orientation', 'horizontal', 'Box', 'off');
    lgd3.Layout.Tile = 'south';

    formatLatex(fig3, "compact-landscape");
    if savePlots
        exportPaperFigure(fig3, ...
            fullfile(outDir, 'field_vs_sim_independent_diagnostics'));
    end
    end

    % Dedicated extended-phi diagnostic, including the proportion of field
    % distances clipped at zero and the rank-matched distance dependence.
    if makeExtendedPhiProfileFigure
    fig4 = figure('Color', 'w', 'Position', [100 100 1000 340]);
    tl4 = tiledlayout(fig4, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax4A = nexttile(tl4, 1);
    yyaxis(ax4A, 'left');
    plot(ax4A, phiProfileValues, distancePhiLoss_m, '-', ...
        'Color', intervalColour, 'LineWidth', 1.8);
    hold(ax4A, 'on');
    plot(ax4A, bestDistancePhi, bestDistanceLoss_m, 'o', ...
        'Color', intervalColour, 'MarkerFaceColor', intervalColour, ...
        'MarkerSize', 6);
    hold(ax4A, 'off');
    ylabel(ax4A, 'Mean quantile gap (m)', 'Interpreter', 'latex');
    xlim(ax4A, [phiProfileValues(1) phiProfileValues(end)]);
    addPhiProfileRangeCue(ax4A, phiOriginalSimulationLimit, ...
        phiProfileValues(end));
    yyaxis(ax4A, 'right');
    plot(ax4A, phiProfileValues, fieldDistanceClippedPct, '--', ...
        'Color', [0.45 0.45 0.45], 'LineWidth', 1.3);
    ylabel(ax4A, 'Zero-clipped field calls (percent)', 'Interpreter', 'latex');
    ylim(ax4A, [0 max(5, 1.08 .* max(fieldDistanceClippedPct))]);
    xlabel(ax4A, 'Effective $\phi$', 'Interpreter', 'latex');
    title(ax4A, '\textbf{A. Extended distance profile}', 'Interpreter', 'latex');
    grid(ax4A, 'on'); grid(ax4A, 'minor'); box(ax4A, 'on');

    ax4B = nexttile(tl4, 2);
    plot(ax4B, phiProfileValues, tbPhiLoss_ms, '-', ...
        'Color', simColour, 'LineWidth', 1.8);
    hold(ax4B, 'on');
    plot(ax4B, bestTbPhi, bestTbLoss_ms, 'o', ...
        'Color', simColour, 'MarkerFaceColor', simColour, ...
        'MarkerSize', 6);
    hold(ax4B, 'off');
    xlabel(ax4B, 'Generative $\phi$', 'Interpreter', 'latex');
    ylabel(ax4B, 'Mean quantile gap (ms)', 'Interpreter', 'latex');
    title(ax4B, '\textbf{B. Extended $T_b$ profile}', 'Interpreter', 'latex');
    xlim(ax4B, [phiProfileValues(1) phiProfileValues(end)]);
    addPhiProfileRangeCue(ax4B, phiOriginalSimulationLimit, ...
        phiProfileValues(end));
    grid(ax4B, 'on'); grid(ax4B, 'minor'); box(ax4B, 'on');

    ax4C = nexttile(tl4, 3);
    fill(ax4C, ...
        [phiDistanceCtr_m(requiredPhiBinKeep); ...
        flipud(phiDistanceCtr_m(requiredPhiBinKeep))], ...
        [requiredPhiQ25(requiredPhiBinKeep); ...
        flipud(requiredPhiQ75(requiredPhiBinKeep))], ...
        fieldColour, 'FaceAlpha', 0.20, 'EdgeColor', 'none');
    hold(ax4C, 'on');
    plot(ax4C, phiDistanceCtr_m(requiredPhiBinKeep), ...
        requiredPhiQ50(requiredPhiBinKeep), '-', ...
        'Color', fieldColour, 'LineWidth', 1.8);
    yline(ax4C, phiOriginalSimulationLimit, 'k--', ...
        'LineWidth', 1.0);
    hold(ax4C, 'off');
    xlabel(ax4C, 'Matched anchor distance (m)', 'Interpreter', 'latex');
    ylabel(ax4C, 'Required effective $\phi$', 'Interpreter', 'latex');
    title(ax4C, '\textbf{C. Distance-binned required $\phi$}', ...
        'Interpreter', 'latex');
    text(ax4C, 0.04, 0.94, ...
        sprintf('global median $=%.2f$', median(requiredPhi, 'omitnan')), ...
        'Units', 'normalized', 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    grid(ax4C, 'on'); grid(ax4C, 'minor'); box(ax4C, 'on');

    formatLatex(fig4, "compact-landscape");
    if savePlots
        exportPaperFigure(fig4, ...
            fullfile(outDir, 'field_vs_sim_phi_profile_extended'));
    end
    end

    if makeD0PhiProfileFigure && runD0PhiSensitivity
        fig5 = figure('Color', 'w', 'Position', [100 100 1000 340]);
        tl5 = tiledlayout(fig5, 1, 3, ...
            'Padding', 'compact', 'TileSpacing', 'compact');

        ax5A = nexttile(tl5, 1);
        imagesc(ax5A, d0PhiProfileValues, d0ProfileValues_m, ...
            d0DistanceLoss_m);
        axis(ax5A, 'xy');
        hold(ax5A, 'on');
        plot(ax5A, d0PhiProfileValues(globalDistancePhiIdx), ...
            d0ProfileValues_m(globalDistanceD0Idx), 'wo', ...
            'MarkerFaceColor', intervalColour, 'MarkerSize', 7, ...
            'LineWidth', 1.0);
        xline(ax5A, phiOriginalSimulationLimit, 'w--', 'LineWidth', 1.0);
        hold(ax5A, 'off');
        xlabel(ax5A, 'Effective $\phi$', 'Interpreter', 'latex');
        ylabel(ax5A, 'Initial distance, $d_0$ (m)', 'Interpreter', 'latex');
        title(ax5A, '\textbf{A. Distance-distribution gap}', ...
            'Interpreter', 'latex');
        cb5A = colorbar(ax5A);
        cb5A.Label.String = 'Mean quantile gap (m)';
        cb5A.Label.Interpreter = 'latex';
        colormap(ax5A, parula);

        ax5B = nexttile(tl5, 2);
        imagesc(ax5B, d0PhiProfileValues, d0ProfileValues_m, ...
            d0TbLoss_ms);
        axis(ax5B, 'xy');
        hold(ax5B, 'on');
        plot(ax5B, d0PhiProfileValues(globalTbPhiIdx), ...
            d0ProfileValues_m(globalTbD0Idx), 'wo', ...
            'MarkerFaceColor', simColour, 'MarkerSize', 7, ...
            'LineWidth', 1.0);
        xline(ax5B, phiOriginalSimulationLimit, 'w--', 'LineWidth', 1.0);
        hold(ax5B, 'off');
        xlabel(ax5B, 'Generative $\phi$', 'Interpreter', 'latex');
        ylabel(ax5B, 'Initial distance, $d_0$ (m)', 'Interpreter', 'latex');
        title(ax5B, '\textbf{B. $T_b$-distribution gap}', ...
            'Interpreter', 'latex');
        cb5B = colorbar(ax5B);
        cb5B.Label.String = 'Mean quantile gap (ms)';
        cb5B.Label.Interpreter = 'latex';
        colormap(ax5B, parula);

        ax5C = nexttile(tl5, 3);
        hDistanceD0 = plot(ax5C, d0ProfileValues_m, ...
            bestDistancePhiByD0, '-o', 'Color', intervalColour, ...
            'MarkerFaceColor', intervalColour, 'LineWidth', 1.8, ...
            'MarkerSize', 5);
        hold(ax5C, 'on');
        hTbD0 = plot(ax5C, d0ProfileValues_m, bestTbPhiByD0, ...
            '-s', 'Color', simColour, 'MarkerFaceColor', simColour, ...
            'LineWidth', 1.8, 'MarkerSize', 5);
        yline(ax5C, phiOriginalSimulationLimit, 'k--', 'LineWidth', 1.0);
        xline(ax5C, defaultSimulationInitialDistance_m, ':', ...
            'Color', [0.35 0.35 0.35], 'LineWidth', 1.0);
        hold(ax5C, 'off');
        xlabel(ax5C, 'Initial distance, $d_0$ (m)', 'Interpreter', 'latex');
        ylabel(ax5C, 'Best-fitting $\phi$', 'Interpreter', 'latex');
        title(ax5C, '\textbf{C. Profile optimum by $d_0$}', ...
            'Interpreter', 'latex');
        xlim(ax5C, [min(d0ProfileValues_m) max(d0ProfileValues_m)]);
        ylim(ax5C, [min(d0PhiProfileValues) max(d0PhiProfileValues)]);
        grid(ax5C, 'on'); grid(ax5C, 'minor'); box(ax5C, 'on');
        lgd5 = legend(ax5C, [hDistanceD0 hTbD0], ...
            {'distance', '$T_b$'}, 'Orientation', 'horizontal', ...
            'Box', 'off');
        lgd5.Location = 'southoutside';

        formatLatex(fig5, "compact-landscape");
        if savePlots
            exportPaperFigure(fig5, ...
                fullfile(outDir, 'field_vs_sim_d0_phi_profile'));
        end
    end

    if makeDistanceBinnedEchoWindowFigure
    fig6 = figure('Color', 'w', 'Position', [100 100 1000 340]);
    tl6 = tiledlayout(fig6, 1, 3, ...
        'Padding', 'compact', 'TileSpacing', 'compact');

    ax6A = nexttile(tl6, 1);
    yyaxis(ax6A, 'left');
    fill(ax6A, ...
        [phiDistanceCtr_m(absoluteTimingBinKeep); ...
        flipud(phiDistanceCtr_m(absoluteTimingBinKeep))], ...
        [requiredPhiQ25(absoluteTimingBinKeep); ...
        flipud(requiredPhiQ75(absoluteTimingBinKeep))], ...
        fieldColour, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    hold(ax6A, 'on');
    plot(ax6A, phiDistanceCtr_m(absoluteTimingBinKeep), ...
        requiredPhiQ50(absoluteTimingBinKeep), '-', ...
        'Color', fieldColour, 'LineWidth', 1.8);
    hold(ax6A, 'off');
    ylabel(ax6A, 'Required effective $\phi$', 'Interpreter', 'latex');
    yyaxis(ax6A, 'right');
    fill(ax6A, ...
        [phiDistanceCtr_m(absoluteTimingBinKeep); ...
        flipud(phiDistanceCtr_m(absoluteTimingBinKeep))], ...
        [callDurationQ25_ms(absoluteTimingBinKeep); ...
        flipud(callDurationQ75_ms(absoluteTimingBinKeep))], ...
        simColour, 'FaceAlpha', 0.16, 'EdgeColor', 'none');
    hold(ax6A, 'on');
    plot(ax6A, phiDistanceCtr_m(absoluteTimingBinKeep), ...
        callDurationQ50_ms(absoluteTimingBinKeep), '-', ...
        'Color', simColour, 'LineWidth', 1.8);
    hold(ax6A, 'off');
    ylabel(ax6A, 'Call duration, $T_c$ (ms)', 'Interpreter', 'latex');
    xlabel(ax6A, 'Matched anchor distance (m)', 'Interpreter', 'latex');
    title(ax6A, '\textbf{A. $\phi$ and call contraction}', ...
        'Interpreter', 'latex');
    set(ax6A, 'XDir', 'reverse');
    grid(ax6A, 'on'); grid(ax6A, 'minor'); box(ax6A, 'on');

    ax6B = nexttile(tl6, 2);
    fill(ax6B, ...
        [phiDistanceCtr_m(absoluteTimingBinKeep); ...
        flipud(phiDistanceCtr_m(absoluteTimingBinKeep))], ...
        [phiWindowQ25_ms(absoluteTimingBinKeep); ...
        flipud(phiWindowQ75_ms(absoluteTimingBinKeep))], ...
        intervalColour, 'FaceAlpha', 0.22, 'EdgeColor', 'none');
    hold(ax6B, 'on');
    plot(ax6B, phiDistanceCtr_m(absoluteTimingBinKeep), ...
        phiWindowQ50_ms(absoluteTimingBinKeep), '-', ...
        'Color', intervalColour, 'LineWidth', 1.8);
    yline(ax6B, 0, 'k:', 'LineWidth', 0.9);
    hold(ax6B, 'off');
    xlabel(ax6B, 'Matched anchor distance (m)', 'Interpreter', 'latex');
    ylabel(ax6B, 'Absolute offset, $\phi T_c$ (ms)', ...
        'Interpreter', 'latex');
    title(ax6B, '\textbf{B. Absolute echo-window term}', ...
        'Interpreter', 'latex');
    set(ax6B, 'XDir', 'reverse');
    grid(ax6B, 'on'); grid(ax6B, 'minor'); box(ax6B, 'on');

    ax6C = nexttile(tl6, 3);
    hDelay6 = plot(ax6C, phiDistanceCtr_m(absoluteTimingBinKeep), ...
        pairedDelayQ50_ms(absoluteTimingBinKeep), '--', ...
        'Color', [0.35 0.35 0.35], 'LineWidth', 1.5);
    hold(ax6C, 'on');
    hPhiWindow6 = plot(ax6C, phiDistanceCtr_m(absoluteTimingBinKeep), ...
        phiWindowQ50_ms(absoluteTimingBinKeep), '-', ...
        'Color', intervalColour, 'LineWidth', 1.8);
    hTa6 = plot(ax6C, phiDistanceCtr_m(absoluteTimingBinKeep), ...
        pairedTaQ50_ms(absoluteTimingBinKeep), '-', ...
        'Color', [0.10 0.10 0.10], 'LineWidth', 1.8);
    hold(ax6C, 'off');
    xlabel(ax6C, 'Matched anchor distance (m)', 'Interpreter', 'latex');
    ylabel(ax6C, 'Duration (ms)', 'Interpreter', 'latex');
    title(ax6C, '\textbf{C. Absolute timing components}', ...
        'Interpreter', 'latex');
    set(ax6C, 'XDir', 'reverse');
    grid(ax6C, 'on'); grid(ax6C, 'minor'); box(ax6C, 'on');
    legend(ax6C, [hDelay6 hPhiWindow6 hTa6], ...
        {'$T_{\mathrm{delay}}$', '$\phi T_c$', '$T_a$'}, ...
        'Location', 'best', 'Box', 'off');

    formatLatex(fig6, "compact-landscape");
    if savePlots
        exportPaperFigure(fig6, ...
            fullfile(outDir, 'field_vs_sim_distance_binned_echo_window'));
    end
    end

    %% Observable agreement and conditional diagnosis
    if makeTimingDiagnosticsFigure
        % Reconstruct the distributions reported numerically at the two
        % conditional profile optima. These summaries remain conditional on
        % fixed k_r and the assumed simulation geometry.
        fieldDistanceAtBestPhi_m = fieldProfileRangeFactor .* max( ...
            fieldProfileTa_s - bestDistancePhi .* fieldProfile.Duration_s, 0);
        simDistanceAtBestPhi_m = SbestPhi.AnchorDistance_m;

        bestTbRows = abs(SphiProfileAll.Phi - bestTbPhi) < 1e-12;
        SbestTbRaw = SphiProfileAll(bestTbRows, :);
        [bestTbSampleRows, bestTbActualCounts] = ...
            sampleRowsByVelocityTarget( ...
            SbestTbRaw.BatSpeed_m_s, profileVelocityMatchEdges_m_s, ...
            profileTargetBinCount, velocityMatchSeed);
        SbestTb = SbestTbRaw(bestTbSampleRows, :);
        simTbAtBestPhi_ms = 1000 .* SbestTb.TbEffective_s;
        fieldTbAtBestPhi_ms = fieldProfileTb_ms;

        diagnosticDistributionSummary = table( ...
            ["Distance"; "Distance"; "T_b"; "T_b"; "Required phi"], ...
            ["Field inferred"; "Simulation"; "Field inferred"; ...
             "Simulation"; "Conditional rank match"], ...
            [bestDistancePhi; bestDistancePhi; bestTbPhi; ...
             bestTbPhi; bestDistancePhi], ...
            [numel(fieldDistanceAtBestPhi_m); ...
             numel(simDistanceAtBestPhi_m); ...
             numel(fieldTbAtBestPhi_ms); numel(simTbAtBestPhi_ms); ...
             numel(requiredPhi)], ...
            [median(fieldDistanceAtBestPhi_m, 'omitnan'); ...
             median(simDistanceAtBestPhi_m, 'omitnan'); ...
             median(fieldTbAtBestPhi_ms, 'omitnan'); ...
             median(simTbAtBestPhi_ms, 'omitnan'); ...
             median(requiredPhi, 'omitnan')], ...
            [prctile(fieldDistanceAtBestPhi_m, 25); ...
             prctile(simDistanceAtBestPhi_m, 25); ...
             prctile(fieldTbAtBestPhi_ms, 25); ...
             prctile(simTbAtBestPhi_ms, 25); ...
             prctile(requiredPhi, 25)], ...
            [prctile(fieldDistanceAtBestPhi_m, 75); ...
             prctile(simDistanceAtBestPhi_m, 75); ...
             prctile(fieldTbAtBestPhi_ms, 75); ...
             prctile(simTbAtBestPhi_ms, 75); ...
             prctile(requiredPhi, 75)], ...
            'VariableNames', {'Metric', 'Dataset', 'ProfilePhi', ...
            'NumCalls', 'Median', 'Q1', 'Q3'});

        diagnosticProfileSummary = table( ...
            ["Distance"; "T_b"], ...
            [bestDistancePhi; bestTbPhi], ...
            [bestDistanceLoss_m; bestTbLoss_ms], ...
            ["m"; "ms"], ...
            [fieldDistanceClippedPct(bestDistanceIdx); NaN], ...
            'VariableNames', {'Profile', 'BestPhi', ...
            'MinimumMeanQuantileGap', 'GapUnit', ...
            'ZeroClippedFieldCalls_pct'});

        fprintf('\nTiming-diagnostic distribution summaries\n');
        disp(diagnosticDistributionSummary);
        disp(diagnosticProfileSummary);
        fprintf('T_b velocity-match shortage across bins: %d calls\n', ...
            sum(profileTargetBinCount - bestTbActualCounts));

        % Create the canvas at its final manuscript dimensions. This keeps
        % the vector bounding box fixed at the same 15.7-cm width as the
        % companion comparison figure, even with an inset or tiled legend.
        figTimingDiagnostics = figure('Color', 'w', 'Units', 'centimeters', ...
            'Position', [2 2 15.70 7.20]);
        tlTimingDiagnostics = tiledlayout(figTimingDiagnostics, 1, 3, ...
            'Padding', 'compact', 'TileSpacing', 'compact');

        % A. Directly observable IPI structure in displacement space.
        axTimingDisplacement = nexttile(tlTimingDiagnostics, 1);
        fieldDispKeep = isfinite(fieldDispCtr) & ...
            isfinite(fieldIpiQ25) & isfinite(fieldIpiQ50) & ...
            isfinite(fieldIpiQ75);
        simDispKeep = isfinite(simDispCtr) & ...
            isfinite(simIpiQ25) & isfinite(simIpiQ50) & ...
            isfinite(simIpiQ75);
        hold(axTimingDisplacement, 'on');
        fill(axTimingDisplacement, ...
            [fieldDispCtr(fieldDispKeep); ...
             flipud(fieldDispCtr(fieldDispKeep))], ...
            [fieldIpiQ25(fieldDispKeep); ...
             flipud(fieldIpiQ75(fieldDispKeep))], ...
            fieldColour, 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        fill(axTimingDisplacement, ...
            [simDispCtr(simDispKeep); flipud(simDispCtr(simDispKeep))], ...
            [simIpiQ25(simDispKeep); flipud(simIpiQ75(simDispKeep))], ...
            simColour, 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        hFieldTiming = plot(axTimingDisplacement, fieldDispCtr(fieldDispKeep), ...
            fieldIpiQ50(fieldDispKeep), '-', 'Color', fieldColour, ...
            'LineWidth', 1.9);
        hSimTiming = plot(axTimingDisplacement, simDispCtr(simDispKeep), ...
            simIpiQ50(simDispKeep), '-', 'Color', simColour, ...
            'LineWidth', 1.9);
        hold(axTimingDisplacement, 'off');
        xlabel(axTimingDisplacement, 'Inter-call displacement (m)', ...
            'Interpreter', 'latex');
        ylabel(axTimingDisplacement, 'IPI (ms)', 'Interpreter', 'latex');
        title(axTimingDisplacement, '\textbf{A. Observable timing}', ...
            'Interpreter', 'latex');
        subtitle(axTimingDisplacement, ...
            {'median and interquartile range', ...
             'binned by inter-call displacement'}, ...
            'Interpreter', 'latex');
        grid(axTimingDisplacement, 'on'); grid(axTimingDisplacement, 'minor'); box(axTimingDisplacement, 'on');
        axis(axTimingDisplacement, 'square');

        % B. Residual timing after exact velocity-distribution matching.
        axTimingVelocity = nexttile(tlTimingDiagnostics, 2);
        fieldVelocityKeep = isfinite(fieldVelocityCtr) & ...
            isfinite(fieldVelocityIpiQ25) & ...
            isfinite(fieldVelocityIpiQ50) & ...
            isfinite(fieldVelocityIpiQ75) & ...
            fieldVelocityBinCount >= diagnosticMinVelocityBinCount;
        simVelocityKeep = isfinite(simVelocityCtr) & ...
            isfinite(simVelocityIpiQ25) & ...
            isfinite(simVelocityIpiQ50) & ...
            isfinite(simVelocityIpiQ75) & ...
            simVelocityBinCount >= diagnosticMinVelocityBinCount;
        hold(axTimingVelocity, 'on');
        fill(axTimingVelocity, ...
            [fieldVelocityCtr(fieldVelocityKeep); ...
             flipud(fieldVelocityCtr(fieldVelocityKeep))], ...
            [fieldVelocityIpiQ25(fieldVelocityKeep); ...
             flipud(fieldVelocityIpiQ75(fieldVelocityKeep))], ...
            fieldColour, 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        fill(axTimingVelocity, ...
            [simVelocityCtr(simVelocityKeep); ...
             flipud(simVelocityCtr(simVelocityKeep))], ...
            [simVelocityIpiQ25(simVelocityKeep); ...
             flipud(simVelocityIpiQ75(simVelocityKeep))], ...
            simColour, 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        plot(axTimingVelocity, fieldVelocityCtr(fieldVelocityKeep), ...
            fieldVelocityIpiQ50(fieldVelocityKeep), '-', ...
            'Color', fieldColour, 'LineWidth', 1.9, ...
            'HandleVisibility', 'off');
        plot(axTimingVelocity, simVelocityCtr(simVelocityKeep), ...
            simVelocityIpiQ50(simVelocityKeep), '-', ...
            'Color', simColour, 'LineWidth', 1.9, ...
            'HandleVisibility', 'off');
        hold(axTimingVelocity, 'off');
        xlabel(axTimingVelocity, 'Velocity (m s$^{-1}$)', ...
            'Interpreter', 'latex');
        ylabel(axTimingVelocity, 'IPI (ms)', 'Interpreter', 'latex');
        title(axTimingVelocity, '\textbf{B. Velocity-conditioned timing}', ...
            'Interpreter', 'latex');
        subtitle(axTimingVelocity, ...
            {'median and interquartile range', ...
             'after velocity-distribution matching'}, ...
            'Interpreter', 'latex');
        xlim(axTimingVelocity, velocityColourLimits_m_s);
        grid(axTimingVelocity, 'on'); grid(axTimingVelocity, 'minor'); box(axTimingVelocity, 'on');
        axis(axTimingVelocity, 'square');

        % C. Normalised conditional distance and T_b loss profiles.
        axTimingProfile = nexttile(tlTimingDiagnostics, 3);
        relativeDistanceLoss = distancePhiLoss_m ./ bestDistanceLoss_m;
        relativeTbLoss = tbPhiLoss_ms ./ bestTbLoss_ms;
        hold(axTimingProfile, 'on');
        hDistanceProfile = plot(axTimingProfile, phiProfileValues, ...
            relativeDistanceLoss, '-', 'Color', intervalColour, ...
            'LineWidth', 1.9);
        hTbProfile = plot(axTimingProfile, phiProfileValues, relativeTbLoss, ...
            '-', 'Color', [0.55 0.34 0.72], 'LineWidth', 1.9);
        plot(axTimingProfile, bestDistancePhi, 1, 'o', ...
            'Color', intervalColour, 'MarkerFaceColor', intervalColour, ...
            'MarkerSize', 5.5, 'HandleVisibility', 'off');
        plot(axTimingProfile, bestTbPhi, 1, 's', ...
            'Color', [0.55 0.34 0.72], ...
            'MarkerFaceColor', [0.55 0.34 0.72], ...
            'MarkerSize', 5.5, 'HandleVisibility', 'off');
        xline(axTimingProfile, phiOriginalSimulationLimit, 'k--', ...
            'LineWidth', 0.9, 'HandleVisibility', 'off');
        hold(axTimingProfile, 'off');
        xlabel(axTimingProfile, 'Echo-window fraction, $\phi$', ...
            'Interpreter', 'latex');
        ylabel(axTimingProfile, 'Relative quantile gap (min. = 1)', ...
            'Interpreter', 'latex'); %keep it min = 1 to fit better
        title(axTimingProfile, '\textbf{C. Acquisition-window diagnostic}', ...
            'Interpreter', 'latex');
        subtitle(axTimingProfile, ...
            {sprintf('dist.: min.-gap $\\phi=%.2f$ ($%.3f$ m)', ...
                bestDistancePhi, bestDistanceLoss_m), ...
             sprintf('$T_b$: min.-gap $\\phi=%.2f$ ($%.2f$ ms)', ...
                bestTbPhi, bestTbLoss_ms)}, ...
            'Interpreter', 'latex');
        xlim(axTimingProfile, [phiProfileValues(1) phiProfileValues(end)]);
        grid(axTimingProfile, 'on'); grid(axTimingProfile, 'minor'); box(axTimingProfile, 'on');
        axis(axTimingProfile, 'square');

        lgdTimingComparison = legend(axTimingDisplacement, [hFieldTiming hSimTiming], ...
            {'field', 'simulation'}, 'Orientation', 'horizontal', ...
            'Box', 'off');
        lgdTimingComparison.Layout.Tile = 'south';
        legend(axTimingProfile, [hDistanceProfile hTbProfile], ...
            {'dist. profile', '$T_b$ profile'}, ...
            'Location', 'northeast', 'Box', 'off');

        % Compact inset confirms the velocity match without using a panel.
        drawnow;
        parentPosition = axTimingVelocity.Position;
        insetPosition = [ ...
            parentPosition(1) + 0.52 * parentPosition(3), ...
            parentPosition(2) + 0.02 * parentPosition(4), ...
            0.44 * parentPosition(3), ...
            0.34 * parentPosition(4)];
        axVelocityInset = axes(figTimingDiagnostics, 'Position', insetPosition, ...
            'Color', 'w', 'Box', 'on');
        fieldVelocitySorted = sort(fieldDiag.Velocity_m_s);
        simVelocitySorted = sort(simDiag.BatSpeed_m_s);
        fieldVelocityCdf = (1:numel(fieldVelocitySorted))' ./ ...
            numel(fieldVelocitySorted);
        simVelocityCdf = (1:numel(simVelocitySorted))' ./ ...
            numel(simVelocitySorted);
        stairs(axVelocityInset, simVelocitySorted, simVelocityCdf, ...
            'Color', simColour, 'LineWidth', 1.0);
        hold(axVelocityInset, 'on');
        stairs(axVelocityInset, fieldVelocitySorted, fieldVelocityCdf, ...
            '--', 'Color', fieldColour, 'LineWidth', 1.0);
        hold(axVelocityInset, 'off');
        xlim(axVelocityInset, velocityColourLimits_m_s);
        ylim(axVelocityInset, [0 1]);
        xticks(axVelocityInset, [0 5 10]);
        yticks(axVelocityInset, [0 1]);
        ylabel(axVelocityInset, 'Cum. frac.', ...
            'Interpreter', 'latex');
        title(axVelocityInset, 'matched velocity', ...
            'Interpreter', 'latex', 'FontWeight', 'normal');
        grid(axVelocityInset, 'on');
        grid(axVelocityInset, 'minor');

        % The compact-wide preset matches the 15.7-cm printed width of the
        % companion figure, preventing LaTeX from shrinking the typography.
        formatLatex(figTimingDiagnostics, "compact-wide");
        axVelocityInset.FontSize = 5.5;
        axVelocityInset.LineWidth = 0.5;
        axVelocityInset.Title.FontSize = 6.0;
        axVelocityInset.XLabel.FontSize = 5.5;
        axVelocityInset.YLabel.FontSize = 5.5;

        if savePlots
            exportPaperFigure(figTimingDiagnostics, ...
                fullfile(outDir, 'field_vs_sim_timing_diagnostics'));
        end
        if saveStats
            writetable(diagnosticDistributionSummary, ...
                fullfile(outDir, ...
                'field_vs_sim_timing_diagnostic_distributions.csv'));
            writetable(diagnosticProfileSummary, ...
                fullfile(outDir, ...
                'field_vs_sim_timing_diagnostic_profiles.csv'));
        end
    end
end

%% Export a cleaned table for downstream modelling
cleanFile = fullfile(outDir, 'vof_processed_data_cleaned.csv');
if saveStats
    writetable(T, cleanFile);
end

fprintf('\nBasic field-data exploration completed.\n');
if savePlots || saveStats || saveSimulationData
    fprintf('Requested outputs saved in: %s\n', outDir);
else
    fprintf('Save switches are off: results were displayed but not written to disk.\n');
end

%% Local helper functions
function flag = robustFlag(x, thresh)
    x = x(:);
    medx = median(x, 'omitnan');
    madx = median(abs(x - medx), 'omitnan');
    if isnan(madx) || madx == 0
        flag = false(size(x));
        return
    end
    robustZ = abs(x - medx) ./ (1.4826 * madx);
    flag = robustZ > thresh;
    flag(isnan(robustZ)) = true;
end

function flag = robustJumpFlag(x, thresh)
    x = x(:);
    dx = [0; diff(x)];
    medx = median(dx, 'omitnan');
    madx = median(abs(dx - medx), 'omitnan');
    if isnan(madx) || madx == 0
        flag = false(size(x));
        return
    end
    robustZ = abs(dx - medx) ./ (1.4826 * madx);
    flag = robustZ > thresh;
    flag(isnan(robustZ)) = true;
    flag(1) = false;
end

function Tcmp = buildTimingComparisonTable(T, S)
    fieldRows = T.IsKeep & isfinite(T.Distance_m) & isfinite(T.Velocity_m_s) & ...
        isfinite(T.IPI_ms) & isfinite(T.Ta_fromIPI_ms) & isfinite(T.Tb_fromIPI_ms);
    simRows = isfinite(S.InterCallDisplacement_m) & isfinite(S.BatSpeed_m_s) & ...
        isfinite(S.IPI_ms) & isfinite(S.Ta_ms) & isfinite(S.Tb_ms);

    Tf = table();
    Tf.Distance_m = T.Distance_m(fieldRows);
    Tf.Velocity_m_s = T.Velocity_m_s(fieldRows);
    Tf.IPI_ms = T.IPI_ms(fieldRows);
    Tf.Ta_ms = T.Ta_fromIPI_ms(fieldRows);
    Tf.Tb_ms = T.Tb_fromIPI_ms(fieldRows);
    Tf.Dataset = categorical(repmat("Field", height(Tf), 1));

    Ts = table();
    Ts.Distance_m = S.InterCallDisplacement_m(simRows);
    Ts.Velocity_m_s = S.BatSpeed_m_s(simRows);
    Ts.IPI_ms = S.IPI_ms(simRows);
    Ts.Ta_ms = S.Ta_ms(simRows);
    Ts.Tb_ms = S.Tb_ms(simRows);
    Ts.Dataset = categorical(repmat("Simulation", height(Ts), 1));

    Tcmp = [Tf; Ts];
end

function outRow = fitFieldSimulationInteractionModel(Tcmp, fieldVar, simVar, responseLabel)
    Tmdl = table();
    Tmdl.Distance_m = Tcmp.Distance_m;
    Tmdl.Velocity_m_s = Tcmp.Velocity_m_s;
    Tmdl.Dataset = reordercats(Tcmp.Dataset, {'Field','Simulation'});
    if strcmp(fieldVar, 'IPI_ms')
        Tmdl.Response_ms = Tcmp.IPI_ms;
    else
        mappedName = matlab.lang.makeValidName(erase(simVar, '_ms')) + "_ms";
        if strcmp(simVar, 'Ta_ms')
            Tmdl.Response_ms = Tcmp.Ta_ms;
        else
            Tmdl.Response_ms = Tcmp.Tb_ms;
        end
    end
    Tmdl = Tmdl(isfinite(Tmdl.Distance_m) & isfinite(Tmdl.Velocity_m_s) & isfinite(Tmdl.Response_ms), :);

    mdl = fitlm(Tmdl, 'Response_ms ~ Distance_m*Velocity_m_s*Dataset');
    coefNames = string(mdl.CoefficientNames(:));
    coefP = mdl.Coefficients.pValue;

    pDataset = pickCoefPValue(coefNames, coefP, "Dataset", false);
    pDistanceByDataset = pickCoefPValue(coefNames, coefP, "Distance_m:Dataset", true);
    pVelocityByDataset = pickCoefPValue(coefNames, coefP, "Velocity_m_s:Dataset", true);
    pThreeWay = pickCoefPValue(coefNames, coefP, "Distance_m:Velocity_m_s:Dataset", true);

    outRow = table( ...
        string(responseLabel), ...
        height(Tmdl), ...
        pDataset, ...
        pDistanceByDataset, ...
        pVelocityByDataset, ...
        pThreeWay, ...
        mdl.Rsquared.Ordinary, ...
        'VariableNames', {'Response','NumRows','P_Dataset','P_DistanceByDataset', ...
        'P_VelocityByDataset','P_DistanceByVelocityByDataset','RSquared'});
end

function pval = pickCoefPValue(coefNames, coefP, pattern, allowContains)
    if allowContains
        idx = contains(coefNames, pattern);
    else
        idx = coefNames == pattern | contains(coefNames, pattern + "_");
    end
    if any(idx)
        pval = max(coefP(idx));
    else
        pval = NaN;
    end
end

function outTbl = summariseResidualSpread(Tcmp, fieldVar, simVar, responseLabel)
    if strcmp(fieldVar, 'IPI_ms')
        responseData = Tcmp.IPI_ms;
    elseif strcmp(simVar, 'Ta_ms')
        responseData = Tcmp.Ta_ms;
    else
        responseData = Tcmp.Tb_ms;
    end

    Tmdl = table(Tcmp.Distance_m, Tcmp.Velocity_m_s, responseData, Tcmp.Dataset, ...
        'VariableNames', {'Distance_m','Velocity_m_s','Response_ms','Dataset'});
    Tmdl.Dataset = reordercats(Tmdl.Dataset, {'Field','Simulation'});
    Tmdl = Tmdl(isfinite(Tmdl.Distance_m) & isfinite(Tmdl.Velocity_m_s) & isfinite(Tmdl.Response_ms), :);
    mdl = fitlm(Tmdl, 'Response_ms ~ Distance_m*Velocity_m_s');
    residuals = mdl.Residuals.Raw;

    outTbl = table();
    for ds = ["Field","Simulation"]
        idx = string(Tmdl.Dataset) == ds;
        thisResidual = residuals(idx);
        thisAbsResidual = abs(thisResidual);
        row = table( ...
            string(responseLabel), string(ds), sum(idx), ...
            median(thisResidual, 'omitnan'), ...
            prctile(thisResidual, 25), prctile(thisResidual, 75), ...
            median(thisAbsResidual, 'omitnan'), ...
            prctile(thisAbsResidual, 75) - prctile(thisAbsResidual, 25), ...
            'VariableNames', {'Response','Dataset','NumRows','MedianResidual_ms', ...
            'ResidualQ1_ms','ResidualQ3_ms','MedianAbsResidual_ms','AbsResidualIQR_ms'});
        outTbl = appendCompatible(outTbl, row);
    end
end

function outTbl = summariseResidualSpreadByDistance(Tcmp, fieldVar, simVar, responseLabel)
    if strcmp(fieldVar, 'IPI_ms')
        responseData = Tcmp.IPI_ms;
    elseif strcmp(simVar, 'Ta_ms')
        responseData = Tcmp.Ta_ms;
    else
        responseData = Tcmp.Tb_ms;
    end

    Tmdl = table(Tcmp.Distance_m, Tcmp.Velocity_m_s, responseData, Tcmp.Dataset, ...
        'VariableNames', {'Distance_m','Velocity_m_s','Response_ms','Dataset'});
    Tmdl.Dataset = reordercats(Tmdl.Dataset, {'Field','Simulation'});
    Tmdl = Tmdl(isfinite(Tmdl.Distance_m) & isfinite(Tmdl.Velocity_m_s) & isfinite(Tmdl.Response_ms), :);
    mdl = fitlm(Tmdl, 'Response_ms ~ Distance_m*Velocity_m_s');
    Tmdl.Residual_ms = mdl.Residuals.Raw;

    edges = linspace(min(Tmdl.Distance_m, [], 'omitnan'), max(Tmdl.Distance_m, [], 'omitnan'), 13);
    binID = discretize(Tmdl.Distance_m, edges);
    centres = 0.5 * (edges(1:end-1) + edges(2:end));

    outTbl = table();
    for ds = ["Field","Simulation"]
        idxDs = string(Tmdl.Dataset) == ds;
        for b = 1:numel(centres)
            idx = idxDs & binID == b;
            if sum(idx) < 5
                continue
            end
            absRes = abs(Tmdl.Residual_ms(idx));
            row = table( ...
                string(responseLabel), string(ds), b, centres(b), sum(idx), ...
                median(absRes, 'omitnan'), prctile(absRes, 25), prctile(absRes, 75), ...
                'VariableNames', {'Response','Dataset','DistanceBin','DistanceMid_m','NumRows', ...
                'MedianAbsResidual_ms','AbsResidualQ1_ms','AbsResidualQ3_ms'});
            outTbl = appendCompatible(outTbl, row);
        end
    end
end

function S = runFieldValidationSimulation(thisDir, krModel, maxVelocityForKeep_m_s, ...
    phiValues, nSequencesPerPhi, initialSpeedRange_m_s, simulationSeedBase, ...
    initialDistance_m)
    addpath(thisDir);

    paramsBase = struct();
    paramsBase.c = 343;
    paramsBase.kr = krModel;
    paramsBase.initialDistance_m = initialDistance_m;
    paramsBase.initialBatSpeed_m_s = 5;
    paramsBase.initialCallDuration_s = 0.010;
    paramsBase.minCallDuration_s = 0.0005;
    paramsBase.maxCalls = 220;
    paramsBase.maxElapsedTime_s = 10;
    paramsBase.maxAnchorDistance_m = 10;
    paramsBase.interceptDistance_m = 0.15;
    paramsBase.numSequences = 1;

    optsBase = struct();
    optsBase.geometryMode = "3D";
    optsBase.numTargets = 3;
    optsBase.targetMotion = true;
    optsBase.targetVelocityMode = "stochastic";
    optsBase.targetVelocityScale = 0.5;
    optsBase.targetVelocityJitterFrac = 0.12;
    optsBase.batVelocityMode = "jittered";
    optsBase.batVelocityJitterFrac = 0.12;
    optsBase.anchorMode = "random";
    optsBase.timingMode = "motionAware";
    optsBase.phiBounds = [0 max(1, max(phiValues))];
    optsBase.callDurationMode = "previousTa";
    optsBase.enforceMaxCallRate = true;
    optsBase.callDurationJitter.enabled = true;
    optsBase.callDurationJitter.mode = "additive";
    optsBase.callDurationJitter.rho = 0.25;

    speedStream = RandStream('mt19937ar', 'Seed', simulationSeedBase);
    seedStream = RandStream('mt19937ar', 'Seed', simulationSeedBase + 1);
    sequenceInitialSpeed_m_s = initialSpeedRange_m_s(1) + ...
        diff(initialSpeedRange_m_s) .* rand(speedStream, nSequencesPerPhi, 1);
    sequenceSimulationSeed = randi(seedStream, 2^31 - 1, ...
        nSequencesPerPhi, 1);

    S = table();
    seqCounter = 0;

    for p = 1:numel(phiValues)
        phi = phiValues(p);
        for s = 1:nSequencesPerPhi
            seqCounter = seqCounter + 1;

            params = paramsBase;
            params.initialBatSpeed_m_s = sequenceInitialSpeed_m_s(s);
            params.maxCallRate_Hz = crMaxFromStop(params, params.initialBatSpeed_m_s);

            opts = optsBase;
            opts.phi = phi;
            opts.rngSeed = sequenceSimulationSeed(s);

            res = simulateResponsivityCore(params, opts);
            Tsim = res.calls;
            if isempty(Tsim)
                continue
            end

            Tsim.SeqID(:) = seqCounter;
            Tsim.Phi = repmat(phi, height(Tsim), 1);
            Tsim.InterCallDisplacement_m = Tsim.BatSpeed_m_s .* Tsim.IPI_s;
            Tsim.IPI_ms = 1000 * Tsim.IPI_s;
            Tsim.Tcall_ms = 1000 * Tsim.Tcall_s;
            Tsim.Ta_ms = 1000 * Tsim.Ta_s;
            Tsim.Tb_ms = 1000 * Tsim.TbEffective_s;

            keepRows = isfinite(Tsim.InterCallDisplacement_m) & isfinite(Tsim.CallRate_Hz) & ...
                isfinite(Tsim.IPI_ms) & isfinite(Tsim.Tcall_ms) & isfinite(Tsim.Ta_ms) & ...
                isfinite(Tsim.Tb_ms) & Tsim.BatSpeed_m_s <= maxVelocityForKeep_m_s;
            Tsim = Tsim(keepRows, :);
            S = appendCompatible(S, Tsim);
        end
    end
end

function [Smatched, summaryTable] = matchSimulationVelocityDistribution( ...
    S, fieldVelocity_m_s, binEdges_m_s, rngSeed)
    % Sample simulated calls within velocity bins to reproduce the retained
    % field call-level velocity distribution without reusing simulations.
    fieldVelocity_m_s = fieldVelocity_m_s(:);
    fieldBin = discretize(fieldVelocity_m_s, binEdges_m_s);
    simBin = discretize(S.BatSpeed_m_s, binEdges_m_s);

    numBins = numel(binEdges_m_s) - 1;
    fieldCount = zeros(numBins, 1);
    availableSimulationCount = zeros(numBins, 1);
    matchedCount = zeros(numBins, 1);
    selectedSimulationRows = false(height(S), 1);

    previousRng = rng;
    cleanupRng = onCleanup(@() rng(previousRng));
    rng(rngSeed);

    for b = 1:numBins
        fieldCount(b) = sum(fieldBin == b);
        candidates = find(simBin == b);
        availableSimulationCount(b) = numel(candidates);
        matchedCount(b) = min(fieldCount(b), availableSimulationCount(b));
        if matchedCount(b) > 0
            chosen = candidates(randperm(numel(candidates), matchedCount(b)));
            selectedSimulationRows(chosen) = true;
        end
    end

    Smatched = S(selectedSimulationRows, :);
    lowerEdge_m_s = binEdges_m_s(1:end-1)';
    upperEdge_m_s = binEdges_m_s(2:end)';
    shortage = fieldCount - matchedCount;
    summaryTable = table(lowerEdge_m_s, upperEdge_m_s, fieldCount, ...
        availableSimulationCount, matchedCount, shortage);

    if any(shortage > 0)
        warning(['Velocity matching could not fill every field bin. ', ...
            '%d field calls lack simulated counterparts.'], sum(shortage));
    end

    clear cleanupRng
end

function [selectedRows, actualBinCount] = sampleRowsByVelocityTarget( ...
    velocity_m_s, binEdges_m_s, targetBinCount, rngSeed)
    % Draw an unrepeated velocity-stratified subset with prescribed counts.
    velocityBin = discretize(velocity_m_s, binEdges_m_s);
    numBins = numel(binEdges_m_s) - 1;
    selectedRows = false(numel(velocity_m_s), 1);
    actualBinCount = zeros(numBins, 1);

    previousRng = rng;
    cleanupRng = onCleanup(@() rng(previousRng));
    rng(rngSeed);

    for b = 1:numBins
        candidates = find(velocityBin == b);
        actualBinCount(b) = min(targetBinCount(b), numel(candidates));
        if actualBinCount(b) > 0
            chosen = candidates(randperm(numel(candidates), actualBinCount(b)));
            selectedRows(chosen) = true;
        end
    end

    clear cleanupRng
end

function [targetBinCount, retainedFieldFraction] = ...
    scaledVelocityTargetCounts(fieldBinCount, availableBinCount)
    % Preserve the field velocity-histogram shape by finding the largest
    % common scalar fraction that every simulation condition can supply.
    fieldBinCount = fieldBinCount(:);
    if size(availableBinCount, 1) ~= numel(fieldBinCount)
        error('Velocity availability rows must match the field bin count.');
    end

    populatedFieldBins = fieldBinCount > 0;
    availabilityRatio = availableBinCount(populatedFieldBins, :) ./ ...
        fieldBinCount(populatedFieldBins);
    retainedFieldFraction = min(1, min(availabilityRatio(:)));
    targetBinCount = floor(retainedFieldFraction .* fieldBinCount);

    minimumAvailability = min(availableBinCount, [], 2);
    restoreOne = populatedFieldBins & targetBinCount == 0 & ...
        minimumAvailability > 0;
    targetBinCount(restoreOne) = 1;
    retainedFieldFraction = sum(targetBinCount) ./ sum(fieldBinCount);
end

function [requiredPhi, pairedAnchorDistance_m, pairedCallDuration_ms, ...
    pairedTa_ms, pairedDelay_ms] = ...
    requiredPhiByVelocityRankMatching(fieldTable, simulationTable, ...
    velocityBinEdges_m_s, krModel, cSound)
    % Within velocity strata, rank field acquisition times against known
    % simulated anchor distances and solve for the effective phi required
    % by each paired observation.
    fieldVelocityBin = discretize(fieldTable.Velocity_m_s, ...
        velocityBinEdges_m_s);
    simVelocityBin = discretize(simulationTable.BatSpeed_m_s, ...
        velocityBinEdges_m_s);
    fieldTa_s = fieldTable.IPI_s ./ (1 + krModel);

    requiredPhi = [];
    pairedAnchorDistance_m = [];
    pairedCallDuration_ms = [];
    pairedTa_ms = [];
    pairedDelay_ms = [];
    for b = 1:(numel(velocityBinEdges_m_s) - 1)
        fieldRows = find(fieldVelocityBin == b);
        simRows = find(simVelocityBin == b);
        numPairs = min(numel(fieldRows), numel(simRows));
        if numPairs == 0
            continue
        end

        [~, fieldOrder] = sort(fieldTa_s(fieldRows), 'ascend');
        [~, simOrder] = sort( ...
            simulationTable.AnchorDistance_m(simRows), 'ascend');
        pairedFieldRows = fieldRows(fieldOrder(1:numPairs));
        pairedSimRows = simRows(simOrder(1:numPairs));

        pairedDistance = simulationTable.AnchorDistance_m(pairedSimRows);
        pairedDelay_s = 2 .* pairedDistance ./ ...
            (cSound + fieldTable.Velocity_m_s(pairedFieldRows));
        phiValues = (fieldTa_s(pairedFieldRows) - pairedDelay_s) ./ ...
            fieldTable.Duration_s(pairedFieldRows);

        requiredPhi = [requiredPhi; phiValues]; %#ok<AGROW>
        pairedAnchorDistance_m = [pairedAnchorDistance_m; pairedDistance]; %#ok<AGROW>
        pairedCallDuration_ms = [pairedCallDuration_ms; ...
            1000 .* fieldTable.Duration_s(pairedFieldRows)]; %#ok<AGROW>
        pairedTa_ms = [pairedTa_ms; ...
            1000 .* fieldTa_s(pairedFieldRows)]; %#ok<AGROW>
        pairedDelay_ms = [pairedDelay_ms; ...
            1000 .* pairedDelay_s]; %#ok<AGROW>
    end
end

function gap = meanAbsoluteQuantileGap(x, y, quantileProbabilities)
    % Robust one-dimensional distributional discrepancy approximating the
    % first Wasserstein distance over the selected quantile interval.
    x = x(isfinite(x));
    y = y(isfinite(y));
    if isempty(x) || isempty(y)
        gap = NaN;
        return
    end

    qPercent = 100 .* quantileProbabilities(:);
    xQuantiles = prctile(x, qPercent);
    yQuantiles = prctile(y, qPercent);
    gap = mean(abs(xQuantiles - yQuantiles), 'omitnan');
end

function names = fieldValidationCacheVariables()
    names = { ...
        'Phi', ...
        'AnchorDistance_m', ...
        'BatSpeed_m_s', ...
        'RelativeClosingVelocityForTiming_m_s', ...
        'CallRate_Hz', ...
        'IPI_s', ...
        'IPI_ms', ...
        'Tcall_s', ...
        'Tcall_ms', ...
        'Ta_ms', ...
        'TbEffective_s', ...
        'Tb_ms', ...
        'InterCallDisplacement_m'};
end

function tableOut = trimFieldValidationSimulationTable(tableIn)
    names = fieldValidationCacheVariables();
    missing = setdiff(names, tableIn.Properties.VariableNames);
    if ~isempty(missing)
        error('exploreFieldData:IncompleteSimulationCache', ...
            'Simulation table is missing required variables: %s', ...
            strjoin(missing, ', '));
    end
    tableOut = tableIn(:, names);
end

function names = durationDistanceCacheVariables()
    names = {'AnchorDistance_m', 'BatSpeed_m_s', 'Tcall_ms'};
end

function tableOut = trimDurationDistanceSimulationTable(tableIn)
    names = durationDistanceCacheVariables();
    missing = setdiff(names, tableIn.Properties.VariableNames);
    if ~isempty(missing)
        error('exploreFieldData:IncompleteDurationCache', ...
            'Duration--distance table is missing required variables: %s', ...
            strjoin(missing, ', '));
    end
    tableOut = tableIn(:, names);
end

function [payload, wasLoaded] = loadSimulationCache( ...
    cacheFile, expectedConfig, loadSavedSimulationData)
    payload = [];
    wasLoaded = false;
    if ~loadSavedSimulationData || ~isfile(cacheFile)
        return
    end

    try
        cached = load(cacheFile, 'payload', 'cacheConfig');
        if isfield(cached, 'payload') && ...
                isfield(cached, 'cacheConfig') && ...
                isequaln(cached.cacheConfig, expectedConfig)
            payload = cached.payload;
            wasLoaded = true;
        else
            fprintf('Ignoring stale simulation cache: %s\n', cacheFile);
        end
    catch cacheError
        warning('exploreFieldData:CacheReadFailed', ...
            'Could not read simulation cache %s (%s). Recomputing.', ...
            cacheFile, cacheError.message);
    end
end

function saveSimulationCache(cacheFile, payload, cacheConfig, ...
    saveSimulationData)
    if ~saveSimulationData
        return
    end

    cacheDir = fileparts(cacheFile);
    if ~exist(cacheDir, 'dir')
        mkdir(cacheDir);
    end
    save(cacheFile, 'payload', 'cacheConfig', '-v7.3');
    fprintf('Saved simulation cache: %s\n', cacheFile);
end

function label = simulationSourceLabel(wasLoaded)
    if wasLoaded
        label = 'saved cache';
    else
        label = 'new deterministic simulation';
    end
end

function result = evaluateDurationEnvelopeCompatibility( ...
    fieldTable, simulationTable, matchSampleCount, samplingSeed, numBins)
    simulationTable = simulationTable( ...
        isfinite(simulationTable.AnchorDistance_m) & ...
        isfinite(simulationTable.Tcall_ms), :);

    if matchSampleCount && height(simulationTable) > height(fieldTable) && ...
            height(fieldTable) > 0
        previousRng = rng;
        cleanupRng = onCleanup(@() rng(previousRng)); %#ok<NASGU>
        rng(samplingSeed);
        keepIdx = randperm(height(simulationTable), height(fieldTable));
        simulationTable = simulationTable(keepIdx, :);
    end

    distanceMaximum_m = max( ...
        [fieldTable.DistanceFar_m; simulationTable.AnchorDistance_m], ...
        [], 'omitnan');
    distanceEdges_m = linspace(0, distanceMaximum_m, numBins);
    [distanceCentres_m, simulationQ05_ms, simulationMedian_ms, ...
        simulationQ95_ms, simulationCount] = binnedQuantiles( ...
        simulationTable.AnchorDistance_m, simulationTable.Tcall_ms, ...
        distanceEdges_m, [5 50 95]);

    distanceCentres_m = distanceCentres_m(:);
    simulationQ05_ms = simulationQ05_ms(:);
    simulationMedian_ms = simulationMedian_ms(:);
    simulationQ95_ms = simulationQ95_ms(:);
    simulationCount = simulationCount(:);

    isCompatible = false(height(fieldTable), 1);
    for i = 1:height(fieldTable)
        overlap = distanceCentres_m >= fieldTable.DistanceNear_m(i) & ...
            distanceCentres_m <= fieldTable.DistanceFar_m(i) & ...
            simulationCount > 0;
        if any(overlap)
            durationLower_ms = min(simulationQ05_ms(overlap), ...
                [], 'omitnan');
            durationUpper_ms = max(simulationQ95_ms(overlap), ...
                [], 'omitnan');
            isCompatible(i) = ...
                fieldTable.Duration_ms(i) >= durationLower_ms && ...
                fieldTable.Duration_ms(i) <= durationUpper_ms;
        end
    end

    result = struct( ...
        'SimulationTable', simulationTable, ...
        'DistanceBinCentre_m', distanceCentres_m, ...
        'SimulationQ05_ms', simulationQ05_ms, ...
        'SimulationMedian_ms', simulationMedian_ms, ...
        'SimulationQ95_ms', simulationQ95_ms, ...
        'SimulationNumCalls', simulationCount, ...
        'IsFieldCompatible', isCompatible, ...
        'CompatibilityPct', 100 * mean(isCompatible, 'omitnan'));
end

function S = runDurationDistanceSimulation(thisDir, krModel, ...
    maxVelocityForKeep_m_s, phiValues, nSequencesPerPhi, simulationSeedBase)
    simulatorFile = fullfile(thisDir, 'simulateResponsivityCore.m');
    if ~isfile(simulatorFile)
        error('Could not find simulateResponsivityCore.m in %s', thisDir);
    end
    addpath(thisDir);

    optsBase = struct();
    optsBase.geometryMode = "3D";
    optsBase.numTargets = 3;
    optsBase.anchorMode = "random";
    optsBase.targetMotion = true;
    optsBase.targetVelocityMode = "stochastic";
    optsBase.targetVelocityScale = 0.5;
    optsBase.targetVelocityJitterFrac = 0.12;
    optsBase.batVelocityMode = "jittered";
    optsBase.batVelocityJitterFrac = 0.12;
    optsBase.timingMode = "motionAware";
    optsBase.callDurationMode = "previousTa";
    optsBase.callDurationJitter.enabled = true;
    optsBase.callDurationJitter.mode = "additive";
    optsBase.callDurationJitter.rho = 0.25;
    optsBase.enforceMaxCallRate = true;

    speedStream = RandStream('mt19937ar', 'Seed', simulationSeedBase);
    seedStream = RandStream('mt19937ar', 'Seed', simulationSeedBase + 1);
    sequenceInitialSpeed_m_s = 4 + ...
        5 .* rand(speedStream, nSequencesPerPhi, 1);
    sequenceSimulationSeed = randi(seedStream, 2^31 - 1, ...
        nSequencesPerPhi, 1);

    Tall = table();
    seqCounter = 0;
    for p = 1:numel(phiValues)
        phi = phiValues(p);
        for s = 1:nSequencesPerPhi
            seqCounter = seqCounter + 1;

            params = struct();
            params.kr = krModel;
            params.initialDistance_m = 5;
            params.stopDistance_m = 0.15;
            params.initialBatSpeed_m_s = sequenceInitialSpeed_m_s(s);
            params.initialCallDuration_s = 0.010;
            params.minCallDuration_s = 0.0005;
            params.maxCalls = 250;
            params.maxElapsedTime_s = 10;
            params.maxAnchorDistance_m = 10;
            params.c = 343;

            if isfield(params, 'stopDistance_m')
                interceptDistance_m = params.stopDistance_m;
            else
                interceptDistance_m = 0.15;
            end
            params.maxCallRate_Hz = (params.c + params.initialBatSpeed_m_s) / ...
                (2 * (1 + params.kr) * interceptDistance_m);

            opts = optsBase;
            opts.phi = phi;
            opts.rngSeed = sequenceSimulationSeed(s);

            res = simulateResponsivityCore(params, opts);
            if isempty(res)
                continue
            end

            if isstruct(res) && isfield(res, 'calls')
                Tsim = res.calls;
            elseif istable(res)
                Tsim = res;
            elseif isstruct(res)
                try
                    Tsim = struct2table(res);
                catch
                    continue
                end
            else
                continue
            end

            if ~istable(Tsim) || isempty(Tsim) || width(Tsim) == 0
                continue
            end

            if ~ismember('Tcall_ms', Tsim.Properties.VariableNames) && ismember('Tcall_s', Tsim.Properties.VariableNames)
                Tsim.Tcall_ms = 1000 * Tsim.Tcall_s;
            end
            if ~ismember('AnchorDistance_m', Tsim.Properties.VariableNames)
                if ismember('AnchorDistance_m', Tsim.Properties.VariableNames)
                    % no-op
                elseif ismember('Target1Distance_m', Tsim.Properties.VariableNames)
                    Tsim.AnchorDistance_m = Tsim.Target1Distance_m;
                elseif ismember('NearestTargetDistance_m', Tsim.Properties.VariableNames)
                    Tsim.AnchorDistance_m = Tsim.NearestTargetDistance_m;
                elseif ismember('InterCallDisplacement_m', Tsim.Properties.VariableNames)
                    Tsim.AnchorDistance_m = Tsim.InterCallDisplacement_m;
                else
                    error('Simulation output does not contain an anchor-distance variable.');
                end
            end

            keepRows = isfinite(Tsim.AnchorDistance_m) & isfinite(Tsim.Tcall_ms) & ...
                isfinite(Tsim.BatSpeed_m_s) & Tsim.BatSpeed_m_s <= maxVelocityForKeep_m_s;
            Tsim = Tsim(keepRows, :);
            if isempty(Tsim)
                continue
            end

            Tsim.SeqID = repmat(seqCounter, height(Tsim), 1);
            Tsim.Phi = repmat(phi, height(Tsim), 1);
            Tall = appendCompatible(Tall, Tsim);
        end
    end
    S = Tall;
end

function [centres, q1, q2, q3, n] = binnedQuantiles(x, y, edges, qVec)
    edges = edges(:);
    centres = edges(1:end-1) + diff(edges) / 2;
    q1 = nan(numel(centres), 1);
    q2 = nan(numel(centres), 1);
    q3 = nan(numel(centres), 1);
    n = zeros(numel(centres), 1);
    for i = 1:numel(centres)
        if i < numel(centres)
            idx = x >= edges(i) & x < edges(i+1);
        else
            idx = x >= edges(i) & x <= edges(i+1);
        end
        vals = y(idx & isfinite(y));
        n(i) = numel(vals);
        if n(i) > 0
            qq = prctile(vals, qVec);
            q1(i) = qq(1);
            q2(i) = qq(2);
            q3(i) = qq(3);
        end
    end
end

function plotDistanceRatePanel(ax, Tsim, velocityCandidates, curveColours, minDistance_m, panelTitle)
    hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on'); grid(ax, 'minor');

    pointVelocityIdx = nan(height(Tsim), 1);
    for ii = 1:height(Tsim)
        [~, pointVelocityIdx(ii)] = min(abs(velocityCandidates - Tsim.BatSpeed_m_s(ii)));
    end

    for ii = 1:numel(velocityCandidates)
        rows = pointVelocityIdx == ii;
        scatter(ax, Tsim.InterCallDisplacement_m(rows), Tsim.CallRate_Hz(rows), 10, ...
            curveColours(ii, :), 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75);
    end

    fitRows = isfinite(Tsim.InterCallDisplacement_m) & isfinite(Tsim.CallRate_Hz) & ...
        Tsim.InterCallDisplacement_m >= minDistance_m & Tsim.CallRate_Hz > 0;
    xFit = linspace(min(Tsim.InterCallDisplacement_m(fitRows)), max(Tsim.InterCallDisplacement_m(fitRows)), 300);
    fitHandles = gobjects(numel(velocityCandidates), 1);
    fitLabels = strings(numel(velocityCandidates), 1);
    for ii = 1:numel(velocityCandidates)
        vRef = velocityCandidates(ii);
        yFit = vRef ./ xFit;
        fitHandles(ii) = plot(ax, xFit, yFit, '-', 'Color', curveColours(ii, :), 'LineWidth', 1.8);
        fitLabels(ii) = sprintf('$v = %.0f$ m s$^{-1}$', vRef);
    end

    xlabel(ax, 'Inter-call displacement (m)', 'Interpreter', 'latex');
    ylabel(ax, 'Call rate (Hz)', 'Interpreter', 'latex');
    title(ax, panelTitle, 'Interpreter', 'latex');
    legend(ax, fitHandles, fitLabels, 'Location', 'northeast');
    hold(ax, 'off');
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

function fitStruct = fitDurationRateHyperbola(callRate_Hz, duration_ms, minCallRate_Hz)
    valid = isfinite(callRate_Hz) & isfinite(duration_ms) & ...
        callRate_Hz >= minCallRate_Hz & duration_ms > 0;

    x = callRate_Hz(valid);
    y = duration_ms(valid);

    fitStruct = struct( ...
        'NumRows', numel(x), ...
        'Offset_ms', nan, ...
        'Scale_msHz', nan, ...
        'RMSE_ms', nan, ...
        'RSquared', nan, ...
        'MinCallRate_Hz', minCallRate_Hz);

    if numel(x) < 3
        return
    end

    X = [ones(numel(x),1), 1 ./ x];
    beta = X \ y;
    yhat = X * beta;

    resid = y - yhat;
    sse = sum(resid .^ 2, 'omitnan');
    sst = sum((y - mean(y, 'omitnan')) .^ 2, 'omitnan');

    fitStruct.Offset_ms = beta(1);
    fitStruct.Scale_msHz = beta(2);
    fitStruct.RMSE_ms = sqrt(mean(resid .^ 2, 'omitnan'));
    fitStruct.RSquared = 1 - sse / sst;
end

function fitStruct = fitDurationRateEnvelope(callRate_Hz, duration_ms, minCallRate_Hz, numBins, q)
    valid = isfinite(callRate_Hz) & isfinite(duration_ms) & ...
        callRate_Hz >= minCallRate_Hz & duration_ms > 0;

    x = callRate_Hz(valid);
    y = duration_ms(valid);

    fitStruct = struct( ...
        'NumRows', numel(x), ...
        'NumBinsUsed', 0, ...
        'Offset_ms', nan, ...
        'Scale_msHz', nan, ...
        'RMSE_ms', nan, ...
        'RSquared', nan, ...
        'MinCallRate_Hz', minCallRate_Hz, ...
        'BinCallRate_Hz', [], ...
        'BinDuration_ms', [], ...
        'Quantile', q);

    if numel(x) < max(10, numBins)
        return
    end

    edges = linspace(min(x), max(x), numBins + 1);
    binX = nan(numBins, 1);
    binY = nan(numBins, 1);

    for i = 1:numBins
        if i < numBins
            inBin = x >= edges(i) & x < edges(i + 1);
        else
            inBin = x >= edges(i) & x <= edges(i + 1);
        end
        if sum(inBin) < 8
            continue
        end
        binX(i) = median(x(inBin), 'omitnan');
        binY(i) = quantile(y(inBin), q);
    end

    keep = isfinite(binX) & isfinite(binY);
    binX = binX(keep);
    binY = binY(keep);
    fitStruct.NumBinsUsed = numel(binX);
    fitStruct.BinCallRate_Hz = binX;
    fitStruct.BinDuration_ms = binY;

    if numel(binX) < 3
        return
    end

    X = [ones(numel(binX),1), 1 ./ binX];
    beta = X \ binY;
    yhat = X * beta;

    resid = binY - yhat;
    sse = sum(resid .^ 2, 'omitnan');
    sst = sum((binY - mean(binY, 'omitnan')) .^ 2, 'omitnan');

    fitStruct.Offset_ms = beta(1);
    fitStruct.Scale_msHz = beta(2);
    fitStruct.RMSE_ms = sqrt(mean(resid .^ 2, 'omitnan'));
    fitStruct.RSquared = 1 - sse / sst;
end

function plotDurationFitCurve(ax, fitStruct, xLimVals, lineSpec, lineWidth, varargin)
    if ~isfinite(fitStruct.Offset_ms) || ~isfinite(fitStruct.Scale_msHz)
        return
    end

    x0 = max([xLimVals(1), fitStruct.MinCallRate_Hz, eps]);
    x1 = xLimVals(2);
    if ~(isfinite(x0) && isfinite(x1) && x1 > x0)
        return
    end

    xFit = linspace(x0, x1, 300);
    yFit = fitStruct.Offset_ms + fitStruct.Scale_msHz ./ xFit;
    if ~isempty(varargin)
        plot(ax, xFit, yFit, lineSpec, 'LineWidth', lineWidth, 'Color', varargin{1});
    else
        plot(ax, xFit, yFit, lineSpec, 'LineWidth', lineWidth);
    end
end

function plotEnvelopeSupportPoints(fitStruct, colourSpec)
    if isempty(fitStruct.BinCallRate_Hz) || isempty(fitStruct.BinDuration_ms)
        return
    end
    plot(fitStruct.BinCallRate_Hz, fitStruct.BinDuration_ms, 'o', ...
        'Color', colourSpec, 'MarkerFaceColor', colourSpec, ...
        'MarkerSize', 4, 'LineWidth', 0.8);
end

function addPhiProfileRangeCue(parentAx, originalLimit, profileMaximum)
    % Shade the phi range used in the original simulations and mark its
    % upper boundary while leaving the extended profile visible.
    yLimits = ylim(parentAx);
    wasHeld = ishold(parentAx);
    hold(parentAx, 'on');
    shade = patch(parentAx, ...
        [0 originalLimit originalLimit 0], ...
        [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
        [0.85 0.88 0.92], 'FaceAlpha', 0.22, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(shade, 'bottom');
    xline(parentAx, originalLimit, 'k--', 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
    xlim(parentAx, [0 profileMaximum]);
    ylim(parentAx, yLimits);
    if ~wasHeld
        hold(parentAx, 'off');
    end
end

function addInlineVelocityScaleStrip(parentAx, valueLimits, xLimVals, yLimVals, labelText)
    hold(parentAx, 'on');

    xSpan = xLimVals(2) - xLimVals(1);
    ySpan = yLimVals(2) - yLimVals(1);

    x0 = xLimVals(1) + 0.63 * xSpan;
    x1 = xLimVals(1) + 0.91 * xSpan;
    y0 = yLimVals(1) + 0.12 * ySpan;
    yLabel = yLimVals(1) + 0.20 * ySpan;
    yTick = yLimVals(1) + 0.07 * ySpan;

    n = 24;
    xVals = linspace(x0, x1, n);
    cVals = linspace(valueLimits(1), valueLimits(2), n);
    scatter(parentAx, xVals, repmat(y0, 1, n), 20, cVals, 's', 'filled', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1.0);

    text(parentAx, x0, yTick, sprintf('%.0f', valueLimits(1)), ...
        'Interpreter', 'latex', 'FontSize', 7, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    text(parentAx, mean([x0 x1]), yTick, sprintf('%.0f', mean(valueLimits)), ...
        'Interpreter', 'latex', 'FontSize', 7, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    text(parentAx, x1, yTick, sprintf('%.0f', valueLimits(2)), ...
        'Interpreter', 'latex', 'FontSize', 7, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    text(parentAx, mean([x0 x1]), yLabel, sprintf('%s (m s$^{-1}$)', labelText), ...
        'Interpreter', 'latex', 'FontSize', 7, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

    hold(parentAx, 'off');
end

function bandIdx = assignNearestBand(values, candidates)
    values = values(:);
    bandIdx = nan(size(values));
    for ii = 1:numel(values)
        [~, bandIdx(ii)] = min(abs(candidates - values(ii)));
    end
end

function [chi2Stat, df, pValue] = chiSquareHomogeneity(obs)
    obs = double(obs);
    rowSums = sum(obs, 2);
    colSums = sum(obs, 1);
    totalN = sum(obs, 'all');
    expected = rowSums * colSums / totalN;
    chi2Stat = sum(((obs - expected) .^ 2) ./ expected, 'all');
    df = (size(obs, 1) - 1) * (size(obs, 2) - 1);
    if exist('chi2cdf', 'file') == 2
        pValue = 1 - chi2cdf(chi2Stat, df);
    else
        pValue = 1 - gammainc(chi2Stat / 2, df / 2);
    end
end
