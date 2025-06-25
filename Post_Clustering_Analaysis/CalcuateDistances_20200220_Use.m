%% 
%Set root and find folders of interest
%root = "Y:\Jui-Yen Huang\Derived Files\Data\In vivo Calcium Imaging\20210129-0208_jRGECO1a_Chrm3q-Ncre19-2";
root = pwd;
%methods = ["", "_Jaccard", "_Cosine"];
methods = [""];
%methods = ["_Cosine"];

for method = methods
    cd(root)
    fullDir_outside = dir;
    dataFoldersLOG_outside = false(length(fullDir_outside),1); % Creates logical array of 0's using 'false'

    for x = 1:length(fullDir_outside)
        dataFoldersLOG_outside(x) = ~isempty(strfind(fullDir_outside(x).name, 'jRGECO1a')) &...
                                    isempty(strfind(fullDir_outside(x).name, 'P10'));
                                    isempty(strfind(fullDir_outside(x).name, 'Summary'));
                             
    end

    folderNames_outside = {fullDir_outside.name}';
    dataFolders_outside = folderNames_outside(dataFoldersLOG_outside);
    dataFolders_outside = dataFolders_outside';
    %%
    summaryPairwiseAll = [];
    summaryCentroidAll = [];
    summaryPWCentroidAll = [];

    summaryPairwiseCluster = [];
    summaryCentroidCluster = [];

    animalNames = {};
    %% Loop through age folders
    for folder_idx = 1:length(dataFolders_outside)
        name = dataFolders_outside{folder_idx};
        cd(dataFolders_outside{folder_idx})
    %%
        fullDir = dir;
        dataFoldersLOG = false(length(fullDir),1); % Creates logical array of 0's using 'false'

        for x = 1:length(fullDir)
               dataFoldersLOG(x) = ~isempty(strfind(fullDir(x).name, 'P')) &...
                               isempty(strfind(fullDir(x).name, 'stim')) &...
                               isempty(strfind(fullDir(x).name, 'Summary'));
            
        end

        folderNames = {fullDir.name}';
        dataFolders = folderNames(dataFoldersLOG);
        dataFolders = dataFolders';
        %% Loop through each animal folder in current age folder
        for x = 1:length(dataFolders)
            try
                load(dataFolders{x} + "/suite2p/plane0/Fall.mat", "iscell", "stat")
            catch ME
                disp("Unable to load " + dataFolders{x} + "/suite2p/plane0/Fall.mat");
                continue
            end
            %% Make the ROI (x,y) locations
            ROI_idx = find(iscell(:,1) == 1);
            ROI_ID = ROI_idx-1;
            stats_for_ROIs = stat(ROI_idx)';

            yx = NaN(length(stats_for_ROIs),2);

            for z = 1:length(yx)
               yx(z,1:2) = stats_for_ROIs{z,1}.med;
            end

            xy = zeros(size(yx, 1), size(yx, 2));
            xy(:,1) = yx(:,2);
            xy(:,2) = yx(:,1);

            xy = round(xy);
            %%
            load_loc = dataFolders{x} + "/Modol_outputs/all" + method + ".mat";
            load_loc = fullfile(load_loc);
            load(load_loc);
            %%
            allCetroidDistances = {}; % Distances of each neuron from its cluster's centroid
            allPairwiseDistances = {}; % Dinstance between each neuron in a cluster
            centroidLocs = [];
            for idx = 1:length(C0)
                centroidDistances = [];
                pairwiseDistances = [];
                clusterLocs = xy(C0{idx}, :);
                centroid = [mean(clusterLocs(:,1)), mean(clusterLocs(:,2))];
                centroidLocs = [centroidLocs; centroid];
                for locIdx = 1:length(clusterLocs)
                    centroidDistances = [centroidDistances; pdist([clusterLocs(locIdx, :); centroid],'euclidean')];
                end

                pairwiseDists = zeros(size(clusterLocs, 1), size(clusterLocs, 1));
                for xIdx = 1:length(clusterLocs)
                    for yIdx = 1:length(clusterLocs)
                        pairwiseDists(xIdx, yIdx) = pdist([clusterLocs(xIdx, :); clusterLocs(yIdx, :)],'euclidean');
                    end
                end
                allCetroidDistances{idx} = centroidDistances;
                allPairwiseDistances{idx} = pairwiseDists;
            end
            % Keeping track of distaces between centroids too to see if there
            % is overlap
            pairwiseCentroidDistances = zeros(size(centroidLocs, 1), size(centroidLocs, 1));
            for xIdx = 1:size(centroidLocs, 1)
                for yIdx = 1:size(centroidLocs, 1)
                    pairwiseCentroidDistances(xIdx, yIdx) = pdist([centroidLocs(xIdx, :); centroidLocs(yIdx, :)],'euclidean');
                end
            end
            Distances.allCentroidDistances = allCetroidDistances;
            Distances.allPairwiseDistances = allPairwiseDistances;
            Distances.pairwiseCentroidDistances = pairwiseCentroidDistances;
            save(dataFolders{x} + "/Modol_outputs/Distances" + method + ".mat", "Distances")

            % All of this was necessary to put the data into Excel files as we
            % discussed. The pattern repeats but essentially the code flattens
            % the data into one column or takes a mean, while also removing
            % duplicates and 0 values
            flatPairwise = [];
            flatCentroid = [];
            tmp = tril(pairwiseCentroidDistances, -1);
            flatPWCentroid = tmp(tmp ~= 0);
            if isempty(C0)
                continue
            end
            for idx = 1:length(C0)
                flatCentroid = [flatCentroid; unique(allCetroidDistances{idx})];     
                tmp = tril(allPairwiseDistances{idx}, -1);
                flatPairwise = [flatPairwise; tmp(tmp ~= 0)];
            end

            x1 = size(summaryPairwiseAll, 1);
            x2 = size(flatPairwise, 1);
            if x2 > x1 && x1 > 0
                summaryPairwiseAll(end+1:end+(x2 - x1), :) = NaN;
            elseif x1 > x2
                flatPairwise(end+1:end+(x1-x2)) = NaN;
            end
            summaryPairwiseAll = [summaryPairwiseAll, flatPairwise];

            x1 = size(summaryCentroidAll, 1);
            x2 = size(flatCentroid, 1);
            if x2 > x1 && x1 > 0
                summaryCentroidAll(end+1:end+(x2 - x1), :) = NaN;
            elseif x1 > x2
                flatCentroid(end+1:end+(x1-x2)) = NaN;
            end
            summaryCentroidAll = [summaryCentroidAll, flatCentroid];

            x1 = size(summaryPWCentroidAll, 1);
            x2 = size(flatPWCentroid, 1);
            if x2 > x1 && x1 > 0
                summaryPWCentroidAll(end+1:end+(x2 - x1), :) = NaN;
            elseif x1 > x2
                flatPWCentroid(end+1:end+(x1-x2), 1) = NaN;
            end
            summaryPWCentroidAll = [summaryPWCentroidAll, flatPWCentroid];


            clusterPairwise = [];
            clusterCentroid = [];
            for idx = 1:length(C0)
                clusterCentroid = [clusterCentroid; mean(allCetroidDistances{idx})];     
                tmp = tril(allPairwiseDistances{idx}, -1);
                clusterPairwise = [clusterPairwise; mean(tmp(tmp ~= 0))];
            end

            x1 = size(summaryCentroidCluster, 1);
            x2 = size(clusterCentroid, 1);
            if x2 > x1 && x1 > 0
                summaryCentroidCluster(end+1:end+(x2 - x1), :) = NaN;
            elseif x1 > x2
                clusterCentroid(end+1:end+(x1-x2), 1) = NaN;
            end
            summaryCentroidCluster = [summaryCentroidCluster, clusterCentroid];

            x1 = size(summaryPairwiseCluster, 1);
            x2 = size(clusterPairwise, 1);
            if x2 > x1 && x1 > 0
                summaryPairwiseCluster(end+1:end+(x2 - x1), :) = NaN;
            elseif x1 > x2
                clusterPairwise(end+1:end+(x1-x2), 1) = NaN;
            end
            summaryPairwiseCluster = [summaryPairwiseCluster, clusterPairwise];

            animalNames{end+1} = dataFolders{x};
        end
        cd("..")
    end
    % Put all the summary matrices together in one cell to loop through
    summaries = {summaryCentroidAll, summaryCentroidCluster, ...
        summaryPWCentroidAll, summaryPairwiseCluster, summaryPairwiseAll};
    filenames = ["summaryCentroidAll", "summaryCentroidCluster", ...
        "summaryPWCentroidAll", "summaryPairwiseCluster", "summaryPairwiseAll"];

    % Loop through and save each summary Excel file
    for idx = 1:length(summaries)
        tab = table([animalNames; num2cell(summaries{idx})]);
        writetable(tab, root+"/Summary_Modol/"+filenames(idx)+method+".xlsx", 'WriteVariableNames', 0, 'WriteMode', 'overwritesheet')
    end
    % Save the .mat data too in case we need it
    save(root+"/Summary_Modol/Distance_summaries" + method + ".mat", "summaries", "filenames","animalNames")
end