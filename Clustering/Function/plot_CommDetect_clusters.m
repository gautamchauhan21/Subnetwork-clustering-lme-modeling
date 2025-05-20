function plot_CommDetect_clusters(all_mat_path, showFigs, saveJPG)
%PLOT_COMMDETECT_CLUSTERS Visualize community detection results with grid-based layout.
%
% Syntax:
%   plot_CommDetect_clusters('path/to/all.mat')
%   plot_CommDetect_clusters('path/to/all.mat', false, true)
%
% Inputs:
%   all_mat_path - Full path to the 'all.mat' file with clustering data
%   showFigs     - (Optional, logical, default = true) Show figures
%   saveJPG      - (Optional, logical, default = true) Save figures as .jpg
%
% Description:
%   This function visualizes community detection results including:
%     - Grid-reordered similarity matrix (based on community index)
%     - Sorted similarity matrix + event raster
%     - Binary raster of significant clusters
%     - Colored raster showing neuron-cluster memberships
%
% Outputs:
%   - Grid_community.fig/.jpg
%   - Event_Clusters.fig/.jpg
%   - Neuronal_Cluster.fig/.jpg
%   - Neuronal_Cluster_colored.fig/.jpg

% == Handle optional inputs
if nargin < 2, showFigs = true; end
if nargin < 3, saveJPG = true; end

% == Output path
[output_path, ~, ~] = fileparts(all_mat_path);

% == Load required data
try
    data = load(all_mat_path);
catch ME
    error('Failed to load data from "%s". Error: %s', all_mat_path, ME.message);
end

% == Figure-specific for community detection: Grid-based communities
% Require function in BCT 
if isfield(data, 'rho') && isfield(data, 'ciu')
    rho   = data.rho;
    ciu = data.ciu;

    fig_grid = figure('Visible', showFigs);
    [gx,gy,idxorder] = grid_communities(ciu); % reorder matrix by community
    imagesc(rho(idxorder,idxorder));
    hold on;
    plot(gx,gy,'k');
    colormap(fcn_cmapjet); % this function is not in BCT so add in the end as helper function
    title('Grid-Reordered Similarity Matrix');
    colorbar;
    hold off;
    
    savefig(fig_grid, fullfile(output_path, 'Grid_community.fig'));
    if saveJPG
        saveas(fig_grid, fullfile(output_path, 'Grid_community.jpg'));
    end
    close(fig_grid);

else
    warning('Skipping Grid_community: Missing rho or ciu.');
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
    warning('Skipping Event_Clusters: Missing MSort, Race, x1, or x2.');
end


% == Figure 2: Raster plot for significant clusters)
if isfield(data, 'RList') && ~isempty(data.RList) && all(isfield(data, {'Race', 'x1'}))
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
     warning('Skipping Neuronal_Cluster: Missing or empty RList, or missing Race/x1.');
end


% == Figure 3: Colored raster plot by neuronal cluster
requiredVars = {'Race', 'RList', 'PCl', 'C0'};
if all(isfield(data, requiredVars)) && ~isempty(data.C0)
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
    warning('Skipping Neuronal_Cluster_colored: Missing Race, RList, PCl, or C0.');
end

end

%% === Helper Function: Custom Colormap
function cmapjet = fcn_cmapjet(len)
    if nargin == 0
        len = 256;
    end
    cmapjet = jet(7);
    cmapjet(4,:) = 1;
    cmapjet = interp1(linspace(0,1,size(cmapjet,1)),cmapjet,linspace(0,1,len));
end
