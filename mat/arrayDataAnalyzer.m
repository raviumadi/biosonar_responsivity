classdef arrayDataAnalyzer
    properties
        AudioData           % Audio data
        SampleRate          % Sampling rate
        FileName            % Input .wav file name
        Segments            % Struct to store segment information
        OutputDir           % Parent directory for storing extracted calls
        MicrophonePositions % Microphone positions
        C                   % Speed of sound
    end

    methods
        % Constructor
        function obj = arrayDataAnalyzer(wavFileName)
            % Load audio data
            [obj.AudioData, obj.SampleRate] = audioread(wavFileName);
            obj.FileName = wavFileName;

            % Create parent directory
            [filePath, fileBase, ~] = fileparts(wavFileName);
            obj.OutputDir = fullfile(filePath, [fileBase '_Analysis']);
            if ~isfolder(obj.OutputDir)
                mkdir(obj.OutputDir);
            end

            % Initialize properties
            obj.Segments = struct('Begin', [], 'End', [], 'Calls', [], 'ValidatedCalls', [], ...
                'Timestamps', [], 'CallRate', [],  'CallDurations', [], ...
                'CallLevels', [], 'BatLocations', [], 'Velocities', []);
            obj.MicrophonePositions = [0 0 0; -1.15 0 -1.1; 0 0 -0.93; 1.15 0 -1.1];
            obj.C = 344; % m/s
            disp('Audio File Loaded Successfully');
        end

        % Select Segments Visually
        function obj = selectSegments(obj)
            figure('Name', 'Select Segments', 'NumberTitle', 'off');
            plot(obj.AudioData);
            title('Click to Select Begin and End of Each Segment');
            xlabel('Sample Number');
            ylabel('Amplitude');

            % Collect segment points
            while true
                [x, ~] = ginput(2); % Two points: Begin and End
                beginPos = round(x(1));
                endPos = round(x(2));

                if beginPos < 1 || endPos > size(obj.AudioData, 1)
                    warning('Selection out of bounds. Try again.');
                    continue;
                end

                segmentIndex = length(obj.Segments) + 1;
                obj.Segments(segmentIndex).Begin = beginPos;
                obj.Segments(segmentIndex).End = endPos;

                % Create subdirectory for the segment
                segmentDir = fullfile(obj.OutputDir, sprintf('segment_%02d', segmentIndex-1));
                if ~isfolder(segmentDir)
                    mkdir(segmentDir);
                end

                % Highlight selected segment
                hold on;
                plot(beginPos:endPos, obj.AudioData(beginPos:endPos), 'r');

                % Ask user if they want to select more segments
                choice = questdlg('Select another segment?', ...
                    'Continue?', ...
                    'Yes', 'No', 'Yes');
                if strcmp(choice, 'No')
                    break;
                end
            end
            close;
        end

        % Set Threshold and Detect Calls for Each Segment
        function obj = setThresholdDetectCalls(obj)
            for i = 2:length(obj.Segments)
                segment = obj.Segments(i);
                segmentData = obj.AudioData(segment.Begin:segment.End, :);
                done = false;
                while ~done
                    % Plot the segment for threshold selection
                    figure;
                    plot(segmentData(:,3));
                    title(sprintf('Set Threshold for Segment %d', i-1));
                    xlabel('Sample Number');
                    ylabel('Amplitude');

                    % Get threshold via mouse click
                    thr = ginput(1);
                    close;
                    threshold = thr(2);

                    % Find calls using the selected threshold
                    calls = obj.findCalls(segmentData(:,3), threshold, 0.005 * obj.SampleRate);

                    % Visualize detected calls
                    figure;
                    plot(segmentData(:,3));
                    hold on;
                    plot(calls, threshold, 'ok');
                    hold off;
                    title(sprintf('Detected Calls for Segment %d', i-1));
                    xlabel('Sample Number');
                    ylabel('Amplitude');

                    % GUI dialog box to confirm threshold
                    choice = questdlg('Is this threshold good?', ...
                        'Threshold Confirmation', ...
                        'Yes', 'No', 'Yes');
                    if strcmp(choice, 'Yes')
                        done = true;
                    end
                    close;
                end

                % Save detected calls for the segment
                obj.Segments(i).Calls = calls;
                segmentDir = fullfile(obj.OutputDir, sprintf('segment_%02d', i-1));
                for j = 1:length(calls)
                    callWindow = max(1, calls(j) - 0.010 * obj.SampleRate): ...
                        min(length(segmentData), calls(j) + 0.010 * obj.SampleRate);
                    callData = segmentData(callWindow,:);
                    filename = sprintf('Call_%02d.wav', j);
                    audiowrite(fullfile(segmentDir, filename), callData, obj.SampleRate);
                end
            end
        end

        % Validate Calls with Spectrogram
        function obj = validateCalls(obj)
            for i = 2:length(obj.Segments)
                calls = obj.Segments(i).Calls;
                validatedCalls = [];
                segmentData = obj.AudioData(obj.Segments(i).Begin:obj.Segments(i).End, :); % Channel 3
                segmentDir = fullfile(obj.OutputDir, sprintf('segment_%02d', i-1));

                for j = 1:length(calls)
                    callIdx = calls(j);
                    callWindow = max(1, callIdx - 0.010 * obj.SampleRate): ...
                        min(length(segmentData), callIdx + 0.010 * obj.SampleRate);
                    callData = segmentData(callWindow,:);

                    % Plot the call and spectrogram
                    figure;
                    subplot(2, 1, 1);
                    plot(callData);
                    title(sprintf('Segment %d - Call %d', i-1, j));
                    xlabel('Sample Number');
                    ylabel('Amplitude');
                    subplot(2, 1, 2);
                    obj.imgSpec(callData(:, 3), min(128, length(callData)), 127, 4096, obj.SampleRate, 100);
                    title('Spectrogram');
                    xlabel('Time (s)');
                    ylabel('Frequency (Hz)');
                    drawnow;

                    % GUI dialog box for user validation
                    choice = questdlg('Keep this call?', ...
                        sprintf('Validation: Segment %d - Call %d', i-1, j), ...
                        'Yes', 'No', 'Yes');
                    if strcmp(choice, 'Yes')
                        validatedCalls = [validatedCalls, calls(j)]; %#ok<AGROW>
                        disp('Call Retained!');
                        % filename = sprintf('Call_%02d.wav', j);
                        % audiowrite(fullfile(segmentDir, filename), callData, obj.SampleRate);
                    else
                        filename = sprintf('Call_%02d.wav', j);
                        delete(fullfile(segmentDir, filename));
                        disp('Call discarded!');
                    end
                    close;
                end

                % Update validated calls in the segment structure
                obj.Segments(i).ValidatedCalls = validatedCalls';
            end
        end

        % Analyze Calls
        function obj = analyzeCalls(obj)
            for i = 2:length(obj.Segments)
                segmentDir = fullfile(obj.OutputDir, sprintf('segment_%02d', i-1));
                callFiles = dir(fullfile(segmentDir, 'Call_*.wav'));

                % Initialize per-segment storage
                obj.Segments(i).CallDurations = [];
                obj.Segments(i).CallLevels = [];
                obj.Segments(i).BatLocations = [];

                for j = 1:length(callFiles)
                    % Read the call
                    callFile = fullfile(callFiles(j).folder, callFiles(j).name);
                    [callData, fs] = audioread(callFile);
                    % callData = obj.cleanCall(callData);

                    % Compute call properties
                    [dur, pres] = obj.percentDurRMS(callData, [0.05, 0.90]);
                    level = 20 * log10(pres / 2e-5);

                    % Cross-correlate to find arrival time
                    delay = zeros(1, size(callData, 2));
                    for k = 1:size(callData, 2)
                        [~, delay(k)] = max(abs(hilbert(xcorr(callData(:, k), callData(:, 2)))));
                    end

                    % Localize bat position
                    batLoc = obj.localization(delay ./ fs, obj.MicrophonePositions, obj.C, '2d');
                    batLoc(:,2) = real(batLoc(:,2)*-1);
                    % batLoc = real(batLoc);
                    % Store in segment structure
                    obj.Segments(i).CallDurations = [obj.Segments(i).CallDurations; dur];
                    obj.Segments(i).CallLevels = [obj.Segments(i).CallLevels; level];
                    obj.Segments(i).BatLocations = [obj.Segments(i).BatLocations; batLoc];
                end
            end
            disp('Call analysis complete.');
        end

        function obj = calculateCallTimestampsVelocityAndRate(obj)
            % Iterate through each segment
            for i = 1:length(obj.Segments)
                segment = obj.Segments(i);
                if isempty(segment.ValidatedCalls) || isempty(segment.BatLocations)
                    continue; % Skip segments with no validated calls or bat locations
                end

                % Extract validated calls and bat locations
                validatedCalls = segment.ValidatedCalls;
                batLocations = segment.BatLocations; % [x, y, z] per call

                % Pre-allocate storage for timestamps and velocities
                numCalls = length(validatedCalls);
                timestamps = zeros(numCalls, 1);
                velocities = zeros(numCalls - 1, 1); % One less velocity than calls

                % Calculate timestamps for each call
                for j = 1:numCalls
                    callSample = validatedCalls(j); % Call sample number
                    timestamps(j) = callSample / obj.SampleRate; % Timestamp in seconds
                end

                % Calculate velocities between consecutive calls
                timeDifferences = diff(timestamps); % Time differences between consecutive calls
                for j = 1:numCalls - 1
                    % Location difference (vector magnitude)
                    locationDiff = norm(batLocations(j + 1, :) - batLocations(j, :));

                    % Velocity (magnitude of movement / time)
                    velocities(j) = locationDiff / timeDifferences(j);
                end

                % Calculate call rate (calls per second)
                if ~isempty(timeDifferences)
                    % avgTimeDiff = mean(timeDifferences); % Average time difference
                    callRate = 1 ./ timeDifferences; % Call rate (calls per second)
                else
                    callRate = 0; % No calls or insufficient data
                end

                % Store results in the segment structure
                obj.Segments(i).Timestamps = timestamps; % Store timestamps
                obj.Segments(i).Velocities = velocities; % Store velocities
                obj.Segments(i).CallRate = callRate;     % Store call rate
            end

            disp('Timestamps, velocities, and call rates calculated and stored in segments.');
        end

        function obj = moveAudioFile(obj)
            % Extract the file path and base name
            [~, fileBase, ~] = fileparts(obj.FileName);

            % Create the destination path in the OutputDir
            destinationPath = fullfile(obj.OutputDir, [fileBase '.wav']);

            % Move the file
            movefile(obj.FileName, destinationPath);

            % Update FileName property to point to the new location. When
            % the object is saved, accessing the file location is possible
            % via the variable.
            obj.FileName = destinationPath;
        end

        % Format and write out the table
        function obj = createAndAppendSegmentsTable(obj, tablePath)
            % 1. Convert Segments struct to table
            segmentsTable = struct2table(obj.Segments);

            % 2. Remove the first row since it's always empty
            if height(segmentsTable) > 1
                segmentsTable(1, :) = []; % Remove first row
            else
                warning('No valid segment data to save. Skipping table update.');
                return;
            end

            % 3. Extract the file base name (without path or extension)
            [~, fileBase, ~] = fileparts(obj.FileName);

            % 4. Add a column with the filename
            segmentsTable.Filename = repmat({fileBase}, height(segmentsTable), 1);

            % 5. Check if the original table exists
            if isfile(tablePath)
                % Load the existing table
                loadedData = load(tablePath, 'AnalysisTable');
                if isfield(loadedData, 'AnalysisTable')
                    existingTable = loadedData.AnalysisTable;

                    % 6. Check if columns match
                    if width(existingTable) ~= width(segmentsTable)
                        error("Column mismatch: Existing table has %d columns, but new data has %d columns.", ...
                            width(existingTable), width(segmentsTable));
                    end

                    % Append new data
                    updatedTable = [existingTable; segmentsTable];

                else
                    warning('MAT file exists but does not contain a valid AnalysisTable. Creating a new table.');
                    updatedTable = segmentsTable;
                end
            else
                % If file doesn't exist, create a new table
                updatedTable = segmentsTable;
            end

            % 7. Save the updated table
            AnalysisTable = updatedTable; %#ok<NASGU>
            save(tablePath, 'AnalysisTable');

            disp(['Segments data successfully appended to: ', tablePath]);
        end
    end

    methods (Access = private)
        % Find calls
        function calls = findCalls(~, signal, threshold, windowLength)
            calls = find(signal > threshold);
            calls = calls([true; diff(calls) > windowLength]);
        end

        % Spectrogram plot
        function imgSpec(~, sig, n, o, F, fs, range)
            [B, f, t] = spectrogram(sig, n, o, F, fs);
            bmin = max(max(abs(B))) / range;
            imagesc(t, f, 20 * log10(max(abs(B), bmin) / bmin));
            set(gca, 'YDir', 'normal');
            colormap(flipud(hot));
        end

        % Percent Duration RMS
        function [dur, pres] = percentDurRMS(~, seg, frac)
            for n = 1:length(seg(1, :))
                [ste(:, n), rd(:, n)] = percentdur(seg(:, n), frac);
            end
            [~, l] = max(rd);
            pres = rms(seg(ste(1, l):ste(2, l), :));
            dur = [ste(1, l), ste(2, l)];
        end

        % Percetn Duration
        function [duration,pres, fmax] = percentdur(seg,frac)
            sm = 2*floor((length(seg)/10)/2)+1;
            csig = cumsum(seg.^2);
            csig(1) = 0;
            csig = csig/max(csig);
            csig = smooth(smooth(csig,sm),sm);
            int = find(csig>frac(1) & csig< frac(2));
            if isempty(int)
                duration = [1 2];
                pres = 0;
            else
                duration = [int(1) int(end)];
                pres = rms(seg(int));
                mag = 20*log10(abs(fft(seg(int(1):int(end)), 192))); % each mag value is 1kHz
                [~, fmax] = max(mag(1:length(mag)/2));

            end
        end

        % Clean Call
        function [call, ipt] = cleanCall(~, sig)
            channels = size(sig, 2);
            call = zeros(size(sig));
            for i = 1:channels
                de = envelope(sig(:, i), 1000, 'rms');
                ipt = findchangepts(de, 'MaxNumChanges', 2, 'Statistic', "rms");
                if length(ipt) == 2
                    sig(ipt(2):end, i) = randn(length(sig(ipt(2):end, i)), 1) / 10000;
                elseif isempty(ipt)
                    call(:, i) = sig(:, i);
                else
                    sig(ipt(1):end, i) = randn(length(sig(ipt(1):end, i)), 1) / 10000;
                end
                call(:, i) = sig(:, i);
            end
        end

        % Localization
        function coords = localization(~, arrivaltime, r, c, dim)
            arrivaltime = arrivaltime(:);
            delay = arrivaltime(2:end) - arrivaltime(1);
            r(1, :) = [];
            if strcmp(dim, '2d')
                A = -2 * [r(:, 1), r(:, 3), c^2 * delay];
                b = ((-r(:, 1).^2) - (r(:, 3).^2) + c^2 * delay.^2);
                m = A \ b;
                coords = [m(1), sqrt((m(3) * c).^2 - m(1).^2 - m(2).^2), m(2)];
            else
                A = -2 * [r(:, 1), r(:, 2), r(:, 3), c^2 * delay];
                b = ((-r(:, 1).^2) - r(:, 2).^2 - (r(:, 3).^2) + c^2 * delay.^2);
                m = A \ b;
                coords = m(1:3)';
            end
        end
    end
end