function summary = reproduce_all(options)
%REPRODUCE_ALL Run the manuscript analyses from one repository entry point.
%
% Examples:
%   reproduce_all
%   reproduce_all(SaveFigures=false, SaveTables=false, RunSSG=false)
%   reproduce_all(KeepFigures=false)
%
% The raw field-recording workflow is interactive and is therefore kept
% separate in field_recording_processing/processFieldRecordings.m. The
% manuscript field analysis begins with data/vof_processed_data.csv.

    arguments
        options.SaveFigures (1,1) logical = true
        options.SaveTables (1,1) logical = true
        options.RunSSG (1,1) logical = true
        options.RunValidation (1,1) logical = true
        options.KeepFigures (1,1) logical = true
        options.StopOnError (1,1) logical = true
    end

    repoRoot = fileparts(mfilename('fullpath'));
    addpath(repoRoot);
    addpath(fullfile(repoRoot, 'fcn'));

    runConfig = struct( ...
        'OverrideOutputSwitches', true, ...
        'SaveFigures', options.SaveFigures, ...
        'SaveTables', options.SaveTables, ...
        'CloseFiguresBeforeRun', ~options.KeepFigures, ...
        'ValidationScenario', "");
    responsivityRunConfig("set", runConfig);
    configCleanup = onCleanup(@() responsivityRunConfig("reset"));

    scripts = [
        "exploreCallRateDistnaceDynamics.m"
        "exploreCallDurationByKr.m"
        "exploreIPIDistancePhiMotile.m"
        "exploreTaDistanceDynamics.m"
        "exploreResponsivityCurve.m"
        "exploreTerminalTbWindow.m"
        "exploreFieldData.m"
    ];
    if options.RunSSG
        scripts(end + 1) = "runCoreSSGAnalysisPipeline.m";
    end

    labels = strings(0, 1);
    status = strings(0, 1);
    elapsed_s = zeros(0, 1);
    messages = strings(0, 1);

    for scriptName = scripts'
        [labels, status, elapsed_s, messages] = runOne( ...
            labels, status, elapsed_s, messages, ...
            repoRoot, scriptName, options.StopOnError);
    end

    if options.RunValidation
        validationScenarios = ["exact", "jitteredBat", "movingTarget"];
        for scenario = validationScenarios
            runConfig.ValidationScenario = scenario;
            responsivityRunConfig("set", runConfig);
            validationFile = fullfile('validation', 'crossCheckCoreCrVr.m');
            label = "validation/crossCheckCoreCrVr.m [" + scenario + "]";
            [labels, status, elapsed_s, messages] = runOne( ...
                labels, status, elapsed_s, messages, ...
                repoRoot, validationFile, options.StopOnError, label);
        end
    end

    summary = table(labels, status, elapsed_s, messages, ...
        'VariableNames', {'Analysis', 'Status', 'Elapsed_s', 'Message'});
    disp(summary);
end

function [labels, status, elapsed_s, messages] = runOne( ...
        labels, status, elapsed_s, messages, repoRoot, scriptName, ...
        stopOnError, displayLabel)

    if nargin < 9
        displayLabel = string(scriptName);
    end
    scriptPath = fullfile(repoRoot, scriptName);
    quotedPath = strrep(scriptPath, '''', '''''');
    fprintf('\n=== %s ===\n', displayLabel);
    timer = tic;
    try
        evalin('base', sprintf('run(''%s'')', quotedPath));
        thisStatus = "completed";
        thisMessage = "";
    catch exception
        thisStatus = "failed";
        thisMessage = string(exception.message);
        warning('reproduce_all:AnalysisFailed', ...
            '%s failed: %s', displayLabel, exception.message);
        if stopOnError
            rethrow(exception);
        end
    end

    labels(end + 1, 1) = displayLabel;
    status(end + 1, 1) = thisStatus;
    elapsed_s(end + 1, 1) = toc(timer);
    messages(end + 1, 1) = thisMessage;
end
