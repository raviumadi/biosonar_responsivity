classdef AudioFileIterator
    properties
        ADS          % AudioDatastore object
        ProcessedLog % File to track processed files
        ProcessedSet % Set of processed file names
    end
    
    methods
        function obj = AudioFileIterator(folderPath, logFile)
            % Constructor: Initialize the audio datastore and load processed files
            if nargin < 2
                logFile = '/Users/ravi/Documents/projects/vo_field/analysis/results/processed_files.txt'; % Default log file
            end
            obj.ProcessedLog = logFile;
            
            % Load already processed files into a set
            obj.ProcessedSet = obj.loadProcessedFiles();
            
            % Create audioDatastore
            obj.ADS = audioDatastore(folderPath, 'IncludeSubfolders', false);
            
             % Filter out files that do not contain "IN" in their names
            obj.ADS.Files = obj.ADS.Files(contains(obj.ADS.Files, 'IN'));

            % Filter out already processed files
            obj.ADS.Files = obj.ADS.Files(~ismember(obj.ADS.Files, obj.ProcessedSet));
        end
        
        function [fileName, obj] = nextFile(obj)
            % Returns the next unprocessed file name in the datastore
            if hasdata(obj.ADS)
                [~, info] = read(obj.ADS);
                fileName = info.FileName;
                
                % Mark file as processed
                obj.markFileAsProcessed(fileName);
            else
                fileName = ''; % No more files
            end
        end
    end
    
    methods (Access = private)
        function processedSet = loadProcessedFiles(obj)
            % Load previously processed files from log
            if exist(obj.ProcessedLog, 'file')
                fid = fopen(obj.ProcessedLog, 'r');
                processedList = textscan(fid, '%s', 'Delimiter', '\n');
                fclose(fid);
                processedSet = string(processedList{1});
            else
                processedSet = string([]);
            end
        end
        
        function markFileAsProcessed(obj, fileName)
            % Append processed file to log
            fid = fopen(obj.ProcessedLog, 'a');
            fprintf(fid, '%s\n', fileName);
            fclose(fid);
        end
    end
end