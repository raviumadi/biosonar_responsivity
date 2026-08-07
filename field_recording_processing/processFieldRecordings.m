function processedFiles = processFieldRecordings(rawDataDir, options)
%PROCESSFIELDRECORDINGS Interactively analyse multichannel field recordings.
%
% This portable entry point consolidates the raw-recording workflow from
% the earlier biosonar_responsivity repository. The original recordings
% are not bundled here. Supply their folder explicitly; all derived calls,
% logs, and MAT tables are written beneath this repository by default.
%
% Example:
%   processFieldRecordings("/path/to/raw/wav", ...
%       Mode="buzzManual", MaxFiles=1)
%
% Modes:
%   "regular"       threshold-assisted call selection
%   "buzzManual"    manually bounded calls (used for terminal sequences)
%   "buzzAutomatic" legacy threshold-assisted buzz workflow

    arguments
        rawDataDir (1,1) string
        options.Mode (1,1) string {mustBeMember(options.Mode, ...
            ["regular", "buzzManual", "buzzAutomatic"])} = "regular"
        options.OutputDir (1,1) string = ""
        options.MaxFiles (1,1) double {mustBePositive} = 1
        options.MoveSourceAudio (1,1) logical = false
        options.FileNameToken (1,1) string = "IN"
    end

    workflowDir = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(workflowDir);
    if options.OutputDir == ""
        outputDir = fullfile(repoRoot, 'results', 'field_recording_processing');
    else
        outputDir = char(options.OutputDir);
    end

    rawDataDir = char(rawDataDir);
    if ~isfolder(rawDataDir)
        error('processFieldRecordings:MissingRawData', ...
            'Raw recording folder does not exist: %s', rawDataDir);
    end
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end

    addpath(workflowDir);
    files = dir(fullfile(rawDataDir, '*.wav'));
    if options.FileNameToken ~= ""
        files = files(contains(string({files.name}), options.FileNameToken));
    end

    logFile = fullfile(outputDir, ...
        sprintf('processed_files_%s.txt', lower(options.Mode)));
    completed = readCompletedLog(logFile);
    fullNames = string(fullfile({files.folder}, {files.name}));
    fullNames = fullNames(~ismember(fullNames, completed));

    if isempty(fullNames)
        fprintf('No unprocessed WAV files matched the selection.\n');
        processedFiles = strings(0, 1);
        return
    end

    nToProcess = min(numel(fullNames), floor(options.MaxFiles));
    processedFiles = strings(0, 1);
    for fileIndex = 1:nToProcess
        wavFile = fullNames(fileIndex);
        fprintf('Analysing %s\n', wavFile);
        [~, fileBase] = fileparts(wavFile);
        extractedDir = fullfile(outputDir, 'extracted_calls', fileBase);

        switch options.Mode
            case "regular"
                analyser = arrayDataAnalyzer(wavFile, extractedDir);
                tableFile = fullfile(outputDir, 'AnalysisTable.mat');
            case "buzzManual"
                analyser = arrayDataAnalyzerwithBuzzManual(wavFile, extractedDir);
                tableFile = fullfile(outputDir, 'AnalysisTable_Buzz.mat');
            case "buzzAutomatic"
                analyser = arrayDataAnalyzerwithBuzz(wavFile, extractedDir);
                tableFile = fullfile(outputDir, 'AnalysisTable_Buzz.mat');
        end

        if ~isfolder(analyser.OutputDir)
            mkdir(analyser.OutputDir);
        end

        analyser = analyser.selectSegments();
        analyser = analyser.setThresholdDetectCalls();
        analyser = analyser.validateCalls();
        analyser = analyser.analyzeCalls();
        analyser = analyser.calculateCallTimestampsVelocityAndRate();
        analyser.createAndAppendSegmentsTable(tableFile);

        if options.MoveSourceAudio
            analyser.moveAudioFile();
        end

        appendCompletedLog(logFile, wavFile);
        processedFiles(end + 1, 1) = wavFile; %#ok<AGROW>
    end
end

function completed = readCompletedLog(logFile)
    if ~isfile(logFile)
        completed = strings(0, 1);
        return
    end
    completed = string(readlines(logFile));
    completed(completed == "") = [];
end

function appendCompletedLog(logFile, wavFile)
    fileId = fopen(logFile, 'a');
    if fileId < 0
        error('processFieldRecordings:LogOpenFailed', ...
            'Could not open processing log: %s', logFile);
    end
    cleanup = onCleanup(@() fclose(fileId));
    fprintf(fileId, '%s\n', wavFile);
end
