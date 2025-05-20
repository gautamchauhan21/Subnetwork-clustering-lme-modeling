function plot_Kmean_clusters(all_mat_path, showFigs, saveJPG)
%PLOT_KMEAN_CLUSTERS Generate clustering figures from all.mat
%
% Syntax:
%   plot_Kmean_clusters('path/to/all.mat')
%   plot_Kmean_clusters('path/to/all.mat', false, true)
%
% Inputs:
%   all_mat_path - (char) Path to the 'all.mat' file produced by run_Kmean
%   showFigs     - (Optional, logical, default = true) Display figures during execution
%   saveJPG      - (Optional, logical, default = true) Save figures as .jpg
%
% Description:
%   This function generates visualizations from K-means clustering output,
%   skipping any figures that cannot be rendered due to missing variables.
%
% Outputs:
%   - Event_Clusters.fig/.jpg (if MSort, Race, x1, x2 exist)
%   - Neuronal_Cluster.fig/.jpg (if Race, RList, x1 exist)
%   - Neuronal_Cluster_colored.fig/.jpg (if Race, RList, PCl, C0 exist)

% == Handle optional arguments
if nargin < 2, showFigs = true; end
if nargin < 3, saveJPG = true; end

% == Get output path
[output_path, ~, ~] = fileparts(all_mat_path);

% == Load required data
try
    data = load(all_mat_path);
catch ME
    error('Failed to load data from "%s". Error: %s', all_mat_path, ME.message);
end

% == Figure 1: Similarity Matrix + Sorted Event Raster
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
     warning('Skipping Event_Clusters figure: missing one or more of [MSort, Race, x1, x2]');
end


% == Figure 2: Raster plot for significant clusters)
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
    warning('Skipping Neuronal_Cluster figure: missing [Race, x1, RList] or RList is empty.');
end


% == Figure 3: Colored raster plot by neuronal cluster
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
    warning('Skipping Neuronal_Cluster_colored figure: missing [Race, RList, PCl, C0] or C0 is empty.');
end

end