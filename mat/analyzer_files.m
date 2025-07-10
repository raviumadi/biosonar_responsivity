%% Generate file list
analysis = "buzz";
% analysis = "regular";

if analysis == "regular"
    iterator = AudioFileIterator('/Users/ravi/Documents/projects/vo_field/data/selected');
    % Define the path for storing the analysis results
    tableFilePath = "/Users/ravi/Documents/projects/vo_field/analysis/results/AnalysisTable.mat";
elseif analysis == "buzz"
    % iterator = AudioFileIterator('/Users/ravi/Documents/projects/vo_field/data/buzz_sequences', '/Users/ravi/Documents/projects/vo_field/analysis/results/processed_files_buzz.txt');
    % tableFilePath = "/Users/ravi/Documents/projects/vo_field/analysis/results/AnalysisTable_Buzz.mat";
    iterator = AudioFileIterator('/Users/ravi/Documents/projects/vo_field/data/responsivity_field_analysis', '/Users/ravi/Documents/projects/vo_field/analysis/results/processed_files_responsivity_analysis.txt');
    tableFilePath = "/Users/ravi/Documents/projects/vo_field/analysis/results/AnalysisTable_Responsivity.mat";
end

%%
continueAnalysis = true;
while continueAnalysis
    % Your file analysis code here
    disp('Analyzing file...');  % Placeholder for actual analysis

    % Get next unprocessed file
    [fileName, iterator] = iterator.nextFile();
    disp(['This file: ', fileName]);

    % Create the analyzer object
    if analysis == "regular"
        analyzer = arrayDataAnalyzer(fileName);
    elseif analysis == "buzz"
        analyzer = arrayDataAnalyzerwithBuzzManual(fileName);
    end

    % Step 1: Select segments
    analyzer = analyzer.selectSegments();

    % Step 2: Set thresholds and detect calls for each segment
    analyzer = analyzer.setThresholdDetectCalls();

    % Step 3: Validate detected calls
    analyzer = analyzer.validateCalls();

    % Step 4: Analyze validated calls
    analyzer = analyzer.analyzeCalls();

    % Step 5: Calculate spatioTemporal Parameters
    analyzer = analyzer.calculateCallTimestampsVelocityAndRate;

    % Step 6: Export results to a table
    analyzer.createAndAppendSegmentsTable(tableFilePath);

    % Step 7: Move the analysed audio file into it's analysis directory
    analyzer.moveAudioFile;
    
    % Ask user if they want to analyze the next file
    choice = questdlg('Analyze next file?', 'Next Iteration', 'Yes', 'No', 'No');
    
    if strcmp(choice, 'No')
        continueAnalysis = false;
    end

end