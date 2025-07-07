%% Batch Clustering Execution Script
% This script loops through session folders, loads detected event data, 
% 
% and applies three clustering methods: k-means, Louvain community detection 
% (uniform), and DBSCAN.
% 
% Clustering is performed using a specified similarity metric (e.g., covariance 
% matrix 'CovM').

%% Initialize workspace
clear
clc
% Add paths to custom function directories for clustering and community detection
% Single quotes (' ') must be used for paths; double quotes (" ") are not supported
addpath 'C:\Users\juiyhuan\OneDrive - Indiana University\Github_Desktop\SubNetwork_search_multilevel2025\Clustering\Function'
addpath 'C:\Users\juiyhuan\OneDrive - Indiana University\Github_Desktop\SubNetwork_search_multilevel2025\Clustering\Function_CommDetec'

%% Initialize parallel computing pool for faster processing
numCores = feature('numcores');
fprintf('This system has %d physical CPU cores available.\n', numCores);

wish_coresToUse = 4;
maxCoresToUse = min(wish_coresToUse, numCores);  
if isempty(gcp('nocreate'))
    fprintf('Starting parallel pool with %d cores...\n', maxCoresToUse);
    parpool; 
else
    currentPool = gcp('nocreate');
    fprintf('Parallel pool already running with %d workers.\n', currentPool.NumWorkers);
end

%% Define root directory and find session folders
root_path = 'C:\Users\juiyhuan\OneDrive - Indiana University\Github_Desktop\SubNetwork_search_multilevel2025\Clustering_Example';

% Identify folders containing relevant session data (e.g., folders with "C2-1" in the name)
folder_list_full = dir(root_path);
folderLog = false(length(folder_list_full),1);
for x = 1: height(folder_list_full)
    folderLog(x) = folder_list_full(x).isdir == 1 &... % check isdir should be folder==1
                   contains(folder_list_full(x).name,'C2-1');
end
folder_list = folder_list_full(folderLog); % Filtered list of session folders

%% Loop over each session folder to run clustering
for run_idx = 1: height(folder_list)
    input_path = fullfile(folder_list(run_idx).folder,folder_list(run_idx).name);
    % Define time per frame (in seconds); used to convert event times to frame indices
    time_per_frame = 0.065;
    
    %% Load input data files
   try
        % Read detected events as a binary matrix (non-zero values converted to 1)
        detected_events = readmatrix(fullfile(input_path,'detected_events.xlsx')) ~= 0;

        % Load start_time vector (try .csv first, fallback to .xlsx)
        try
            start_time = readmatrix(fullfile(input_path,'start_time.csv'));
        catch
            start_time = readmatrix(fullfile(input_path,'start_time.xlsx'));
        end
        
        % Load end_time vector (try .csv first, fallback to .xlsx)
        try
            end_time = readmatrix(fullfile(input_path,'end_time.csv'));
        catch
            end_time = readmatrix(fullfile(input_path,'end_time.xlsx'));
        end

    catch
        % Display error if input files are missing or unreadable
        warning('Required input files (start_time, end_time, detected_events) are missing or unreadable.');
        continue;
    end
    
    %% Extract time windowed data for clustering
    % Convert start/end times (in seconds) to frame indices
    % Crop the detected_events matrix accordingly
    Race = detected_events(:, floor(start_time/time_per_frame):floor(end_time/time_per_frame)-1);

    fprintf("process %s \n", folder_list(run_idx).name);
    
    tic  % Start timing this session's clustering run
    run_Kmean(input_path, Race, 'CovM', 100, 5000);
    run_CommDetect_uniform(input_path, Race, 'CovM', 500, 5000);
    run_DBSCAN(input_path, Race,'CovM',5000);
   
    toc  % End timing
    
end