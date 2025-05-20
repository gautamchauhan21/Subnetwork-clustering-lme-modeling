function plot_DBSCAN_clusters(all_mat_path, showFigs, saveJPG)
%PLOT_DBSCAN_CLUSTERS Visualize DBSCAN clustering results from all.mat output.
%
% Syntax:
%   plot_DBSCAN_clusters('path/to/all.mat')
%   plot_DBSCAN_clusters('path/to/all.mat', false, true)
%
% Inputs:
%   all_mat_path - (char) Full path to the 'all.mat' file from run_DBSCAN
%   showFigs     - (Optional, logical, default = true) Show figures during execution
%   saveJPG      - (Optional, logical, default = true) Save figures as .jpg images
%
% Description:
%   This function visualizes DBSCAN clustering outputs including:
%     - k-distance graph for epsilon estimation
%     - Cluster label scatter plot
%     - Sorted similarity matrix and event raster
%     - Binary raster for significant clusters
%     - Colored neuron-cluster raster
%
% Outputs:
%   - k_Distance_graph.fig/.jpg
%   - DBSCAN_clustering_results.fig/.jpg
%   - Event_Clusters.fig/.jpg
%   - Neuronal_Cluster.fig/.jpg
%   - Neuronal_Cluster_colored.fig/.jpg

% == Handle optional arguments
if nargin < 2, showFigs = true; end
if nargin < 3, saveJPG = true; end

% == Output directory
[output_path, ~, ~] = fileparts(all_mat_path);

% == Load data
try
    data = load(all_mat_path);
catch ME
    error('Failed to load data from "%s". Error: %s', all_mat_path, ME.message);
end

% == Figure: k-distance graph
if all(isfield(data, {'sortedDist', 'kneeIdx', 'epsilon', 'k'}))
    k     = data.k;
    sortedDist = data.sortedDist;
    kneeIdx    = data.kneeIdx;
    epsilon    = data.epsilon;
    fig_kDistance = figure('Visible', showFigs);
    hold on
    plot(sortedDist, 'b-', 'LineWidth', 1.5);
    plot(kneeIdx, epsilon, 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Highlight the knee point
    xlabel('Points sorted by distance');
    ylabel(sprintf('%d-th Nearest Neighbor Distance',k));
    title('k-Distance Graph');
    grid on;
    hold off
    
    savefig(fig_kDistance, fullfile(output_path, 'k_Distance_graph.fig'));
    if saveJPG
       saveas(fig_kDistance, fullfile(output_path, 'k_Distance_graph.jpg'));
    end
    close(fig_kDistance);
else
    warning('Skipping k-Distance graph: missing sortedDist, kneeIdx, epsilon, or k.');

end

% == Figure: DBSCAN Clustering Label Scatter Plot
if isfield(data, 'IDX2')
    IDX2  = data.IDX2;

    fig_DBSCAN = figure('Visible', showFigs);
    scatter(1:length(IDX2), IDX2, 15, IDX2, 'filled');
    colormap('jet');
    xlabel('Frame Index'); 
    ylabel('Cluster Label');
    title('DBSCAN Clustering Results');
    grid on;
    
    savefig(fig_DBSCAN, fullfile(output_path, 'DBSCAN_clustering_results.fig'));
    if saveJPG
        saveas(fig_DBSCAN, fullfile(output_path, 'DBSCAN_clustering_results.jpg'));
    end
    close(fig_DBSCAN);
else
    warning('Skipping DBSCAN clustering figure: IDX2 not found.');
end

% == Figure Similarity Matrix + Sorted Event Raster
if all(isfield(data, {'MSort', 'Race', 'x1', 'x2'}))
    
    Race  = data.Race;
    MSort = data.MSort;
    x1    = data.x1;
    x2    = data.x2;
    
    fig1 = figure('Visible', showFigs);
    subplot(2,1,1)
    imagesc(MSort)
    colormap jet
    axis image
    xlabel('RACE #'), ylabel('RACE #')
    title('Sorted Similarity Matrix')
    colorbar
    
    subplot(2,1,2)
    imagesc(Race(x1,x2), [-1 1.2])
    xlabel('RACE #'), ylabel('Neuron #')
    title(' Raster Plot- Sorted Event by cluster');
    
    savefig(fig1, fullfile(output_path, 'Event_Clusters.fig'));
    if saveJPG
        saveas(fig1, fullfile(output_path, 'Event_Clusters.jpg'));
    end
    close(fig1);
else
    warning('Skipping Event_Clusters figure: missing one or more of MSort, Race, x1, x2.');
end


% == Figure Raster plot for significant clusters)
if all(isfield(data, {'Race', 'x1', 'RList'})) && ~isempty(data.RList)
    RList = data.RList;
    
    fig2 = figure('Visible', showFigs);
    imagesc(Race(x1, RList) ~= 0)
    colormap(flipud(gray))
    title('Sorted Raster plot: significant cluster only)')
    
    savefig(fig2, fullfile(output_path, 'Neuronal_Cluster.fig'));
    if saveJPG
        saveas(fig2, fullfile(output_path, 'Neuronal_Cluster.jpg'));
    end
    close(fig2);

else
    warning('Skipping Neuronal_Cluster figure: missing Race, x1, RList or RList is empty.');
end

% == Figure Colored raster plot by neuronal cluster
if all(isfield(data, {'Race', 'RList', 'PCl', 'C0'})) && ~isempty(data.C0)
    C0    = data.C0;
    PCl   = data.PCl;

    % Define color map
    newCmap = [255 255 255; 227 26 28; 31 120 180; 178 223 138; 51 160 44;
               251 154 153; 166 206 227; 253 191 111; 255 127 0; 202 178 214;
               106 61 154; 255 255 153; 177 89 40; 0 0 0;];
    newCmap = newCmap ./ 255;
    
    % Build sorted neuron list
    sortInds = [];
    for ndx = 1:length(C0)
        sortInds = [sortInds, C0{ndx}];
    end
    
    for ndx = 1:size(Race, 1)
        if ~ismember(ndx, sortInds)
            sortInds = [sortInds, ndx];
        end
    end
            
    Race_sort = double(Race(sortInds,RList));
    
    PCl_sort = PCl(:, RList);
    
    % Assign cluster number on rastermap (so colors show up in imagesc)
    for jdx = 1:size(PCl_sort, 2)
        for kdx = 1:length(C0)
            neurs = C0{kdx};
            sortedNeurs = [];
            for neur = neurs
               sortedNeurs = [sortedNeurs, find(sortInds == neur)];
            end
            if PCl_sort(kdx, jdx) == 1
                for mdx = sortedNeurs
                    if Race_sort(mdx, jdx)
                        Race_sort(mdx, jdx) = kdx + 1;
                    end
                end
            end
        end
    end

    fig3 = figure('Visible', showFigs);
    imagesc(Race_sort)
    scMap = [0 0 0; newCmap(1:max(Race_sort(:)), :)];
    
    scMap(2, :) = 0.3;
    
    colormap(scMap)
    
    xlabel("Sorted GCE")
    ylabel("Contour #")
    
    title('Raster Plot: Sorted event and neurons')
    
    savefig(fig3, fullfile(output_path, 'Neuronal_Cluster_colored.fig'));
    if saveJPG
        saveas(fig3, fullfile(output_path, 'Neuronal_Cluster_colored.jpg'));
    end
    close(fig3);

else
    warning('Skipping Neuronal_Cluster_colored figure: missing Race, RList, PCl, C0 or C0 is empty.');
end

end