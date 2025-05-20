function run_DBSCAN(input_path, Race, method,NShuff_num)
%RUN_DBSCAN Perform DBSCAN clustering on event-based neural data
%
% Syntax:
%   run_DBSCAN(input_path, Race, method, NShuff_num)
%
% Inputs:
%   input_path   - (char) Path to session folder (used for saving outputs)
%   Race         - [NCell x NFrame] binary matrix of neural events (1 = event)
%   method       - (char) Similarity method: 'CovM', 'CosineM', or 'JaccardM'
%   NShuff_num   - (int) Number of shuffles for neuron significance testing (e.g., 5000)
%
% Description:
%   This function performs clustering on neural activity using DBSCAN. It estimates
%   the epsilon parameter automatically using a k-distance graph, applies DBSCAN
%   with a minimum point threshold of 3, and processes the results to identify
%   statistically significant neuron-cluster relationships.
%
%   Cluster participation across time and neuron membership are assessed using
%   shuffled null distributions. The final results are saved for downstream
%   visualization and analysis.
%
% Output Files:
%   Results are saved to:
%       [input_path]/Output_DBSCAN_<method>/
%
%   Saved output files include:
%     - all.mat           : All variables
%     - output_result.mat : Summary statistics in struct `output_result`
%
% Output (in `output_result` struct):
%   output_result.total_cell_number         : Total number of neurons
%   output_result.NCl_beforeStat            : Number of event clusters before statistical filtering
%   output_result.silhs_mean_beforeStat     : Mean silhouette score of event clusters
%   output_result.NCl                       : Final number of statistically significant neuron ensembles
%   output_result.No_assemblies (%)         : Percent of timepoints with no ensembles
%   output_result.S_assemblies (%)          : Percent of timepoints with one ensemble
%   output_result.M_assemblies (%)          : Percent of timepoints with more than one ensemble
%   output_result.Cells_not_in_assembly (%) : Percent of neurons in no ensemble
%   output_result.Cells_in_one_assembly (%) : Percent of neurons in one ensemble
%   output_result.Cells_in_many_assembly (%) : Percent of neurons in >1 ensemble
%
% Example:
%   run_DBSCAN('session01', RaceMatrix, 'CosineM', 5000);


%% Compute similarity matrix
switch method
    case 'CovM'
        CovM = CovarM(Race);
    case 'CosineM'
        CovM = CosineM(Race);
    case 'JaccardM'
        CovM = JaccardM(Race);
    otherwise
        error('Unknown method "%s". Valid options: CovM, CosineM, Jaccard.', method);
end

%% Create output folder
output_path = fullfile(input_path, "Output_DBSCAN_" + method);
if ~exist(output_path, 'dir')
    mkdir(output_path)
end

%% Convert Correlation to distance matrix
distMatrix = 1 - CovM; % Convert correlation to distance
distMatrix = distMatrix / max(distMatrix(:)); % Normalize distances to [0, 1]

%%  K-Distance graph for Epsilon estimation
k = 5; % Choose k
kDist = zeros(size(distMatrix, 1), 1);
for i = 1:size(distMatrix, 1)
    sortedD = sort(distMatrix(i, :));
    kDist(i) = sortedD(k); % Distance to the k-th nearest neighbor
end
sortedDist = sort(kDist);

% Knee point detection
startPoint = [1, sortedDist(1)];
endPoint = [length(sortedDist), sortedDist(end)];
distances = zeros(size(sortedDist));
for i = 1:length(sortedDist)
    point = [i, sortedDist(i)];
    distances(i) = abs(det([endPoint - startPoint; point - startPoint])) / norm(endPoint - startPoint);
end

% Find optimal epsilon below threshold
remainingPeaks = distances; % Copy distances to track remaining peaks
validEpsilon = false;
%lastCheckedEpsilon = NaN; % Initialize to store the last checked epsilon

while ~validEpsilon
    [~, kneeIdx] = max(remainingPeaks); % Find the current maximum peak
    epsilon = sortedDist(kneeIdx);     % Corresponding epsilon value
    lastCheckedEpsilon = epsilon;      % Store the last checked epsilon

    if epsilon <= 0.25
        validEpsilon = true; % Stop if epsilon is valid
    else
        remainingPeaks(kneeIdx) = -inf; % Remove this peak and find the next one
        if all(remainingPeaks == -inf)
            disp('No valid epsilon found below 0.25.'); % Fail-safe if no valid epsilon exists
            epsilon = 0.1; % Use the last checked epsilon
            validEpsilon = true; % Exit the loop after setting epsilon
        end
    end
end


%% Apply DBSCAN by using the epsilon identifed with K distance graph
minPts = 3; % using the minimum points
labels = dbscan(distMatrix, epsilon, minPts, 'Distance', 'precomputed');

if all(isnan(labels), 'all')
    warning("DBSCAN could not assign any cluster. Saving empty results.");
    save(fullfile(output_path, 'Clusters.mat'), '-v7.3');
    return
end

% Analyze Clusters as needed
numClusters = numel(unique(labels(labels > 0))); % Exclude noise
%disp(['Number of clusters: ', num2str(numClusters)]);
%disp(['Noise points: ', num2str(sum(labels == -1))]);

clustered_idx = find(labels ==1);
if isempty(clustered_idx)
    disp("no cluster")
    C0 = {}; % need this when summarizing data from all session
    save(fullfile(output_path, 'all.mat'), '-v7.3');
    return
end

%% Below start using the output from DBSCAN and analyzed by the rest of Model code
% first need to convert the DBSCAN output simialr to Kmean output
[NCell,NRace] = size(Race); %total cell/neuron number; total frame number

IDX2 = labels; % Cluster labels from DBSCAN
NCl = max(IDX2, [], 'all'); % Determine the number of clusters (including noise)

% Handle noise points (-1): Assign them to a new cluster index (NCl + 1)
noise_idx = find(IDX2 == -1);
if ~isempty(noise_idx)
    IDX2(IDX2 == -1) = NCl + 1; % Reassign noise points to a new cluster
    NCl = NCl + 1; % Update the total number of clusters
end

% Silhouette values 
silhV_all = silh(CovM, IDX2); % Silhouette values for all clusters

% Calculate median silhouette values for each cluster
sCl_all = zeros(1, NCl); % Preallocate array for silhouette values
for i = 1:NCl
    sCl_all(i) = median(silhV_all(IDX2 == i));
end

% Remove the silhouette value of the noise cluster (last cluster)
if ~isempty(noise_idx)
    sCl_all = sCl_all(1:end-1);
end

% Sort silhouette values in descending order
[sCl_all, xCl] = sort(sCl_all, 'descend');
sCl = rmmissing(sCl_all); % Remove NaN values, if any --silv of event cluster

% Determine the number of valid clusters after sorting
NCl_beforeStat = length(sCl); % Event cluster number before assigning cell

% === Format similar to Kmeans
[~,x2] = sort(IDX2);
MSort = CovM(x2,x2); % Sorted normalized covariance matrix

% === Cluster scores
R = cell(0);
CellScore = zeros(NCell,NCl);
CellScoreN = zeros(NCell,NCl);
for i = 1:NCl
    R{i} = find(IDX2==i);
    CellScore(:,i) = sum(Race(:,R{i}),2);
    CellScoreN(:,i) = CellScore(:,i)/length(R{i});
end

%Assign cells to cluster with which it most likely spikes
[~,CellCl] = max(CellScoreN,[],2);
%Remove cells with less than 2 spikes in a given cluster
CellCl(max(CellScore,[],2)<2) = 0;
[X1,x1] = sort(CellCl);


%% Neuron-to-Cluster Significance Testing
NShuf = NShuff_num; % lower for testing; use 5000 for full

% Count number of participation to each cluster
CellP = zeros(NCell,NCl);
CellR = zeros(NCell,NCl);

parfor i = 1:NCl
    CellP(:,i) = sum(Race(:,IDX2 == i),2);
    CellR(:,i) = CellP(:,i)/sum(IDX2 == i);
end

% Test for statistical significance
CellCl = zeros(NCl,NCell); % Binary matrix of cell associated to clusters
for j = 1:NCell
    %Random distribution among Clusters
    RClr = zeros(NCl,NShuf);
    Nrnd = sum(Race(j,:) ~= 0);
    if Nrnd == 0
        continue
    end
    
    for l = 1:NShuf
        Random = randperm(NRace);
        Random = Random(1:Nrnd);
        Racer = zeros(1,NRace);
        Racer(Random) = 1;
        for i = 1:NCl
            RClr(i,l) = sum(Racer(:,IDX2 == i),2);
        end
    end
    RClr = sort(RClr,2);
  
    %Proba above 95th percentile
    ThMax = RClr(:,floor(NShuf*(1-0.05/NCl))); 
    for i = 1:NCl
        CellCl(i,j) = double(CellP(j,i)>ThMax(i));% - double(RCl(:,j)<ThMin);
    end
end

A0 = find(sum(CellCl) == 0); %Cells not in any cluster
A1 = find(sum(CellCl) == 1); %Cells in one cluster
A2 = find(sum(CellCl) >= 2); %Cells in several clusters

% Resolve multi-cluster assignments
for i = A2
    [~,idx] = max(CellR(i,:));
    CellCl(:,i) = 0;
    CellCl(idx,i) = 1;
end

C0 = cell(0);    % idx of neurons in each subnetwork
k = 0;
inds = [];
for i = 1:NCl
    if length(find(CellCl(i,:)))>2
        k = k+1;
        C0{k} = find(CellCl(i,:));
        inds = [inds; k];
    end
end

% Participation rate to its own cluster
CellParticip = max(CellR([A1 A2],:),[],2);

%% Cluster participation across time
NCl = length(C0); %final subnetwork number
if ~NCl
    NCl = 0;
    disp('There were no significant clusters found!!');
    save(fullfile(output_path, 'all.mat'), '-v7.3');    
    return
    %exit() % in hpc job need use exit    
end

%Cell count in each cluster
RCl = zeros(NCl,NRace);
PCl = zeros(NCl,NRace);
for i = 1:NCl
    RCl(i,:) = sum(Race(C0{i},:));
end

RCln = zeros(NCl,NRace);
for j = 1:NRace
    %Random distribution among Clusters
    RClr = zeros(NCl,NShuf);
    Nrnd = sum(Race(:,j) ~= 0); 
    if ~Nrnd % neuron doesn't fire in time period
        continue
    end

    for l = 1:NShuf
        Random = randperm(NCell);
        Random = Random(1:floor(Nrnd));
        Racer = zeros(NCell,1);
        Racer(Random) = 1;
        for i = 1:NCl
            RClr(i,l) = sum(Racer(C0{i}));
        end
    end
    
    RClr = sort(RClr,2);
    
    %Proba above 95th percentile
    ThMax = RClr(:,round(NShuf*(1-0.05/NCl)));
    
    for i = 1:NCl
        PCl(i,j) = double(RCl(i,j)>ThMax(i));
    end
    %Normalize (probability)
    RCln(:,j) = RCl(:,j)/sum(Race(:,j));
end

%% Final preparations for downstream plotting
if ~NCl
    disp('There were no significant clusters found!!');
end

% Times that significantly recruit 0 cell assemblies; will not plot this
Cl0 = find(sum(PCl,1) == 0);
% Times that significantly recruit 1 cell assembly
Cl1 = find(sum(PCl,1) == 1);
% etc.
Cl2 = find(sum(PCl,1) == 2);
Cl3 = find(sum(PCl,1) == 3);
Cl4 = find(sum(PCl,1) == 4);

Bin = 2.^(0:NCl-1);

%Sort Cl1
[~,x01] = sort(Bin*PCl(:,Cl1));
Cl1 = Cl1(x01);

%Sort Cl2
[~,x02] = sort(Bin*PCl(:,Cl2));
Cl2 = Cl2(x02);

%Sort Cl3
[~,x03] = sort(Bin*PCl(:,Cl3));
Cl3 = Cl3(x03);

RList = [Cl1 Cl2 Cl3 Cl4];
[~,x1] = sort(Bin*CellCl(inds, :));

% Save everything for plotting and result
save(fullfile(output_path, 'all.mat'), '-v7.3');

%% Generate the output_result structure
output_result.total_cell_number = NCell; % total neuronal number
output_result.NCl_beforeStat = length (sCl); % number of event cluster
output_result.silhs_mean_beforeStat = mean(sCl);% mean of silhouette value of event cluster
output_result.NCl = NCl;% number of subnetwork/ensemble
% Total number of GCE times
total = size ([Cl0, Cl1, Cl2, Cl3, Cl4], 2);

% Event Times that recruit no assemblies
if total > 0
    output_result.No_assemblies = length(Cl0)/total*100;
else
    output_result.No_assemblies = NaN; 
end

% Event Times that recruit one assembly
if total > 0
    output_result.S_assemblies = length(Cl1)/total*100;
else
    output_result.S_assemblies = NaN; 
end

% Event Times that recruit more than one assembly
if total > 0
    output_result.M_assemblies = length([Cl2, Cl3, Cl4])/total*100;
else
    output_result.M_assemblies = NaN;
end

% number of neurons that participate in one, more than one or no assemblies
if NCell > 0
    output_result.Cells_not_in_assembly = length(A0)/NCell*100;
    output_result.Cells_in_one_assembly = length(A1)/NCell*100;
    output_result.Cells_in_many_assembly =length(A2)/NCell*100;
else
    output_result.Cells_not_in_assembly = NaN;
    output_result.Cells_in_one_assembly = NaN; 
    output_result.Cells_in_many_assembly = NaN;

% Save output_data for latter summary
save(fullfile(output_path, 'output_result.mat'),'output_result');

end
