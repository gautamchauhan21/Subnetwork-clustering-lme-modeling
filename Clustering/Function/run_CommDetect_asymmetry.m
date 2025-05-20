function run_CommDetect_asymmetry(input_path, Race, method, nreps_num, NShuff_num)
%RUN_COMMDETECT_ASYMMETRY Perform asymmetric community detection on neural data
%
% Syntax:
%   run_CommDetect_asymmetry(input_path, Race, method, nreps_num, NShuff_num)
%
% Inputs:
%   input_path   - (char) Path to session folder (used for saving outputs)
%   Race         - [NCell x NFrame] binary matrix of neural events (1 = event)
%   method       - (char) Similarity method: 'CovM', 'CosineM', or 'JaccardM'
%   nreps_num    - (int) Number of repetitions for modularity optimization per gamma
%                   (e.g., 500 for CosineM, 100 for CovM)
%   NShuff_num   - (int) Number of shuffles for neuron significance testing (e.g., 5000)
%
% Description:
%   This function performs community detection using the asymmetric null model
%   described in Rubinov & Sporns (2011). It handles positive and negative correlations
%   separately using a "negative_asym" model during modularity optimization. It runs
%   multiple gamma values, computes consensus clustering, and derives final community
%   assignments using `consensus_und`.
%
%   Neuron-to-cluster assignments are validated statistically via shuffle-based
%   null models. The function outputs key matrices and summary metrics required for
%   visualization and quantitative comparisons.
%
% Output Files:
%   Results are saved to:
%       [input_path]/Output_CommDetect_Asymmetry_<method>/
%
%   Saved output files include:
%     - all.mat           : All variables
%     - output_result.mat : Summary statistics in struct `output_result`
%
% Output (in `output_result` struct):
%   output_result.total_cell_number          : Total number of neurons
%   output_result.NCl_beforeStat             : Number of event clusters before statistical filtering
%   output_result.silhs_mean_beforeStat      : Mean silhouette score of event clusters
%   output_result.NCl                        : Final number of statistically significant neuron ensembles
%   output_result.No_assemblies (%)          : Percent of timepoints with no ensemble recruitment
%   output_result.S_assemblies (%)           : Percent of timepoints recruiting one ensemble
%   output_result.M_assemblies (%)           : Percent of timepoints recruiting multiple ensembles
%   output_result.Cells_not_in_assembly (%)  : Percent of neurons not in any ensemble
%   output_result.Cells_in_one_assembly (%)  : Percent of neurons in one ensemble
%   output_result.Cells_in_many_assembly (%) : Percent of neurons in multiple ensembles
%
% Example:
%   run_CommDetect_asymmetry('session01', RaceMatrix, 'CosineM', 500, 5000);

%% Compute similarity matrix
switch method
    case 'CovM'
        rho = CovarM(Race);
    case 'CosineM'
        rho = CosineM(Race);
    case 'JaccardM'
        rho = JaccardM(Race);
    otherwise
        error('Unknown method "%s". Valid options: CovM, CosineM, Jaccard.', method);
end

%% Create output folder
output_path = fullfile(input_path, "Output_CommDetect_Asymmetry_" + method);
if ~exist(output_path, 'dir')
    mkdir(output_path)
end


%gamma is a parameter that one can vary--for community size
%In our case, it is difficult to decide which gamma should, so we run different gamma value (gammavals)
gammavals = linspace(0,1,11); % use different gamma 

nreps = nreps_num;
%Pre-allocate matrix / array for data storage 
% store the ci from each given gamma value (each vector) and repeat (colume in vector) into the array.
ci_all = cell(1,length(gammavals)); 
% store the AgreementMatrix from each given gamma (each vector) into the array.
Agreement_all = cell(1,length(gammavals));
% store the q from each given gamma value (row) and repeat (colume) into the array.
q_all = zeros(length(gammavals),nreps);
% store the agreement matrix from all given gamma value and their repeat.
d = zeros(length(rho));
% store the agreement matrix from a given gamma value for its repeat, temp for each given gamma value.
Agreement_temp = zeros(length(rho));

for igamma = 1:length(gammavals) % use different gamma for community_louvain
    gamma = gammavals(igamma);
     
    ci = zeros(length(rho),nreps);
    q = zeros(1,nreps);
    parfor irep = 1:nreps
           [ci(:,irep),q(:,irep)] = community_louvain(rho,[],[],'negative_asym');
    end
    
    % calculate the agreement matrix for the repeat of given gamma value
    Agreement_temp = Agreement_temp + agreement(ci); 
    
    % generate agreement matrix for all repeat from all given gamma value.
    d = d + agreement(ci);
    
    %save data from the given gamma value into array/vector            
    ci_all{igamma} = ci;
    q_all(igamma,:) = q;
    Agreement_all{igamma} = Agreement_temp;
                      
end

%% using the concensus matrix
%normalize the concensus matrix into probability
d_sum = sum (d,"all");
d_Pro = d/d_sum;
d_Pro_max = max (d_Pro,[], 'all');
% The value of d_Pro is small, thus using tau = 0, not setting threshold
% Using the ciu from consensus_und function to assign community
ciu = consensus_und(d_Pro,0,10);

%% Below start using the output from community detection and analyzed by the rest of Model code
%assign the variable for further step

[NCell,NRace] = size(Race); %total cell/neuron number; total frame number
IDX2 = ciu;
IDX2 = transpose (IDX2);

% communit/cluster number
NCl= max(IDX2); % this is the initial community / cluster number
%%
% to make the code identical to Kmean, assign CovM to rho
CovM = rho;
%Calculate silhValue of each cluster/community after community detection before statistic test
silhV_all = silh(CovM, IDX2); % Silhouette values for all clusters

sCl_all = zeros(1,NCl); % to store the silV of each cluster
for i = 1:NCl
    sCl_all(i) = median(silhV_all(IDX2==i));
end

[sCl_all,xCl] = sort(sCl_all,'descend'); %sort the silhValue of all cluster/community
sCl = rmmissing (sCl_all); % remove NaN; silv of event cluster
NCl_beforeStat = length(sCl); % Event cluster number before assigning cell
%%
%Taking the index output from community detection and use that follow the rest of Modol code
[~,x2] = sort(IDX2);

% Sorted normalized covariance matrix

MSort = CovM(x2,x2);

% detected events clusters and their scores
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

for i = 1:NCl
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

C0 = cell(0);   % idx of neurons in each subnetwork
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
    disp('There were no significant clusters found!! Cannot run this cell...');
    save(fullfile(output_path, 'all.mat'),'-v7.3');    
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
    disp('There were no significant clusters found!! Cannot run this cell...');
    save(fullfile(output_path, 'all.mat'), '-v7.3');
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
