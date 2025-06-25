% Clear the workspace and command window
clear
clc
%%
% Add paths to custom function directories for clustering and community detection
% Single quotes (' ') must be used for paths; double quotes (" ") are not supported
addpath 'C:\Users\juiyhuan\OneDrive - Indiana University\Github_Desktop\SubNetwork_search_multilevel2025\Clustering\Function'
addpath 'C:\Users\juiyhuan\OneDrive - Indiana University\Github_Desktop\SubNetwork_search_multilevel2025\Clustering\Function_CommDetec'
%%
% Initialize parallel computing pool for faster processing
if isempty(gcp('nocreate'))
    disp('Starting parallel pool...');
    parpool;  % Start parallel pool with default profile and number of workers
else
    disp('Parallel pool already running.');
end
%%
% Define the input folder path containing session data
input_path = 'C:\Users\juiyhuan\OneDrive - Indiana University\Github_Desktop\SubNetwork_search_multilevel2025\Clustering_Example\P15_C2-1_BR_L1_150';

% Define time per frame (in seconds) for converting time to frame indices
time_per_frame = 0.065;

%% Load input data
% Attempt to read input files (detected_events, start_time, end_time) from the input path
try
    % Read detected events as a binary matrix (non-zero values converted to 1)
    detected_events = readmatrix(fullfile(input_path,'detected_events.xlsx')) ~= 0;
    % Get start_time
    start_time = readmatrix(fullfile(input_path,'start_time.csv'));
    % Get end_time
    end_time = readmatrix(fullfile(input_path,'end_time.csv'));
catch
    % Display error if input files are missing or unreadable
    error('Required input files (start_time, end_time, detected_events) are missing or unreadable.');
end

%%
%% Prepare data for clustering
% Convert time (sec) to frame index using the sampling rate (time_per_frame)
% Crop the event matrix to the time window of interest
Race = detected_events(:, floor(start_time/time_per_frame):floor(end_time/time_per_frame)-1);
%%
%% ------------------------ K-means Clustering ------------------------
% Run K-means clustering using 3 different similarity metrics
% Inputs:
%   input_path   - (char) Path to session folder (used for saving outputs)
%   Race         - [NCell x NFrame] binary matrix of neural events (1 = event)
%   method       - (char) Similarity method: 'CovM', 'CosineM', or 'JaccardM'
%   iterationN   - (int) Number of iterations for K-means optimization
%   NShuff_num   - (int) Number of shuffles for neuron significance testing (e.g., 5000)

run_Kmean(input_path, Race, 'CovM', 10, 100);
run_Kmean(input_path, Race, 'CosineM', 10, 100);
run_Kmean(input_path, Race, 'JaccardM', 10, 100);

% Plot results from Covariance-based K-means clustering
% Define path to the output file from K-means clustering (CovM metric)
all_mat_path = fullfile(input_path,"Output_Kmean_CovM/all.mat");
plot_Kmean_clusters(all_mat_path);
%%
%% ------------------ Community Detection Clustering ------------------
% Run Louvain community detection in both uniform and asymmetric resolution modes
% Inputs:
%   input_path   - (char) Path to session folder (used for saving outputs)
%   Race         - [NCell x NFrame] binary matrix of neural events (1 = event)
%   method       - (char) Similarity method: 'CovM', 'CosineM', or 'JaccardM'
%   nreps_num    - (int) Number of repetitions for modularity optimization per gamma
%                   (e.g., 100 for CosineM, 500 for CovM)
%   NShuff_num   - (int) Number of shuffles for neuron significance testing (e.g., 5000)

run_CommDetect_uniform(input_path, Race, 'CovM', 100, 100);
run_CommDetect_asymmetry(input_path, Race, 'CovM', 100, 100);

% Define path to the output file from uniform Community Detection (CovM metric)
all_mat_path = fullfile(input_path,"Output_CommDetect_Uniform_CovM/all.mat");

% Plot the Community Detection clustering results using the output file
plot_CommDetect_clusters(all_mat_path);

% Define path to the output file from asymmetrical Community Detection (CovM metric)
all_mat_path = fullfile(input_path,"Output_CommDetect_Asymmetry_CovM/all.mat");
% Plot the Community Detection clustering results using the output file
plot_CommDetect_clusters(all_mat_path);
%%
%% ------------------------ DBSCAN Clustering ------------------------
% Run DBSCAN clustering with covariance-based similarity
% Inputs:
%   input_path   - (char) Path to session folder (used for saving outputs)
%   Race         - [NCell x NFrame] binary matrix of neural events (1 = event)
%   method       - (char) Similarity method: 'CovM', 'CosineM', or 'JaccardM'
%   NShuff_num   - (int) Number of shuffles for neuron significance testing (e.g., 5000)

run_DBSCAN(input_path, Race,'CovM',1000);
% Define path to the output file from DBSCAN clustering (CovM metric)
all_mat_path = fullfile(input_path,"Output_DBSCAN_CovM/all.mat");

% Plot the DBSCAN clustering results using the output file
plot_DBSCAN_clusters(all_mat_path);
%% 
%