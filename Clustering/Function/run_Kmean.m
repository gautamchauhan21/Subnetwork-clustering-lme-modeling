function run_Kmean(input_path, Race, method, iterationN, NShuff_num, runNull)
%RUN_KMEAN Perform K-means clustering on event-based neural data.
%
% Syntax:
%   run_Kmean(input_path, Race, method, iterationN, NShuff_num)
%   run_Kmean(input_path, Race, method, iterationN, NShuff_num, runNull)
%
% Inputs:
%   input_path   - (char) Path to session folder (used for saving outputs)
%   Race         - [NCell x NFrame] binary matrix of neural events (1 = event)
%   method       - (char) Similarity method: 'CovM', 'CosineM', or 'JaccardM'
%   iterationN   - (int) Number of iterations for K-means optimization
%   NShuff_num   - (int) Number of shuffles for neuron significance testing (e.g., 5000)
%   runNull      - (Optional, logical, default = false) Run null model clustering
%
% Description:
%   This function performs unsupervised K-means clustering on neural event data
%   based on similarity of activation patterns. It supports optional null model testing
%   to validate cluster significance. The output includes cluster labels, sorted
%   similarity matrices, neuron-to-cluster assignments, and ensemble participation data.
%
%   The function saves results in a subfolder of `input_path` named:
%     "Output_Kmean_<method>"
%
%   Saved output files include:
%     - all.mat           : All variables
%     - output_result.mat : Summary statistics in struct `output_result`
%
% Output (in output_result structure):
%   output_result.total_cell_number        : Total number of neurons
%   output_result.NCl_beforeStat           : Number of clusters before stat filtering
%   output_result.silhs_mean_beforeStat    : Mean silhouette score before stat
%   output_result.NCl                      : Final number of statistically significant clusters
%   output_result.No_assemblies (%)        : Percent of timepoints with no ensembles
%   output_result.S_assemblies (%)         : Percent of timepoints with one ensemble
%   output_result.M_assemblies (%)         : Percent of timepoints with multiple ensembles
%   output_result.Cells_not_in_assembly (%) : Percent of neurons in no ensemble
%   output_result.Cells_in_one_assembly (%) : Percent of neurons in one ensemble
%   output_result.Cells_in_many_assembly (%) : Percent of neurons in >1 ensemble
%
% Example:
%   run_Kmean('session01', RaceMatrix, 'CovM', 100, 5000);
%   run_Kmean('session01', RaceMatrix, 'CosineM', 100, 5000, true);  % With null model test


if nargin < 6
    runNull = false;
end

%% Init & Size
[NCell, NRace] = size(Race);

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
output_path = fullfile(input_path, "Output_Kmean_" + method);
if ~exist(output_path, 'dir')
    mkdir(output_path)
end

%% Run k-means clustering
[IDX2, sCl, ~, ~] = kmeansopt(CovM, iterationN, 'precomp');
NCl = max(IDX2);
[~,x2] = sort(IDX2);
MSort = CovM(x2,x2); % Sorted normalized covariance matrix

%% Run null model clustering (random trials)(optional)
if runNull
    fprintf('Running null model clustering...\n')
    NTrials = 1000; % Use 1000 for sufficient repeat
    
    sClrnd = zeros(1,NTrials);
    for i = 1:NTrials
        sClrnd(i) = kmeansoptrnd(Race,100,NCl,'var');
    end
    
    NClOK = sum(sCl>prctile(sClrnd(1:NTrials-1),95));
    sClOK = sCl(1:NClOK)';
    
    RaceOK = Race(:,IDX2<=NClOK);
    NRaceOK = size(RaceOK,2);
    
    save(fullfile(output_path, 'Null_model_Clusters.mat'));
end

%% Score calculation for each cell and cluster
R = cell(1, NCl);
CellScore = zeros(NCell,NCl);
CellScoreN = zeros(NCell,NCl);
for i = 1:NCl
    R{i} = find(IDX2==i);
    CellScore(:,i) = sum(Race(:,R{i}),2);
    CellScoreN(:,i) = CellScore(:,i)/length(R{i});
end

%% Assign cells to cluster with which it most likely spikes
[~,CellCl] = max(CellScoreN,[],2);
CellCl(max(CellScore,[],2)<2) = 0; % remove weakly spiking cells

[X1,x1] = sort(CellCl); % for later plotting

NCl_beforeStat = length(sCl); % Event cluster number before assigning cell


%% Neuron-to-Cluster Significance Testing
NShuf = NShuff_num; % lower for testing; use 5000 for full

% Count number of participation to each cluster
CellP = zeros(NCell,NCl);
CellR = zeros(NCell,NCl);

for i = 1:NCl
    CellP(:,i) = sum(Race(:,IDX2 == i),2);
    CellR(:,i) = CellP(:,i)/sum(IDX2 == i);
end

% Test for statistical significance
CellCl = zeros(NCl,NCell); % Binary matrix of cell associated to clusters
for j = 1:NCell
    % Random distribution among Clusters
    RClr = zeros(NCl,NShuf);
    Nrnd = sum(Race(j,:) ~= 0);
    if Nrnd == 0
        continue
    end
    
    parfor l = 1:NShuf
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
NCl = length(C0);  %final subnetwork number
if ~NCl
    NCl = 0;
    warning('No significant cluster found in %s',input_path);
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
    % Random distribution among Clusters
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
    warning('No significant clusters found in %s',input_path);
    return;
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
end

% Save output_data for latter summary
save(fullfile(output_path, 'output_result.mat'),'output_result');

end
