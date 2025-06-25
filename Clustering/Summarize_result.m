%%  Summary of clustering analysis results from multiple sessions


% Clear the workspace and command window
clear
clc
%%
%% Define path and search for session folders
% Set root directory containing clustering result folders
root_path = 'C:\Users\juiyhuan\OneDrive - Indiana University\Github_Desktop\SubNetwork_search_multilevel2025\Clustering_Example';

% Identify folders containing relevant session data (e.g., folders with "C2-1" in the name)
folder_list_full = dir(root_path);
folderLog = false(length(folder_list_full),1);
for x = 1: height(folder_list_full)
    folderLog(x) = folder_list_full(x).isdir == 1 &... % check isdir should be folder==1
                   contains(folder_list_full(x).name,'C2-1');
end
folder_list = folder_list_full(folderLog);

%% Preallocate variables to store metadata and result summaries
% Metadata
all_name = strings(0,1);
all_LitterID = strings(0,1);
all_age = [];
all_SubjectID = strings(0,1);
all_location = [];
all_Layer = [];

% Clustering output summary
all_total_cell_number = [];
all_NCl_beforeStat = [];
all_silhs_mean_beforeStat = [];
all_NCl = [];
all_No_assemblies = [];
all_S_assemblies = [];
all_M_assemblies = [];
all_Cells_not_in_assembly = [];
all_Cells_in_one_assembly = [];
all_Cells_in_many_assembly = [];

%% Loop over each session folder to extract clustering result and metadata
for idx = 1: height(folder_list)
    % Construct full path to current session
    session_path = fullfile(folder_list(idx).folder,folder_list(idx).name);
    
    % Define expected output result file path
    result_path = fullfile(session_path, 'Output_Kmean_CovM','output_result.mat');
    
    % Skip folder if result file is missing
    if ~exist(result_path,'file')
        warning("No output_result.mat file in %s",session_path);
    else            
        try
            load(result_path);  % Load 'output_result' structure
        catch ME
            warning ("Unable to load result at %s",result_path); continue;
        end
    end
    
    %% Extract metadata from folder name
    session_name = string(folder_list(idx).name);
    meta_infor = split(session_name,'_');
    age_infor = char(meta_infor{1,1});
    LitterID_infor = string(meta_infor{2,1});
    subject_infor = string(meta_infor{3,1});
    image_location_infor = string(meta_infor{4,1});
    cortical_layer_infor = str2double(meta_infor{end,1});

    age = str2double(age_infor(2:end)); % Remove 'P' prefix to get numeric age
    Litter_ID = strcat("C57-",LitterID_infor);
    subject_id = strcat("C57-",LitterID_infor, "_", subject_infor);
    
    % Identify cortical layer from z-position threshold based on age
    if age == 11 || age == 13 || age == 15 
        
        if cortical_layer_infor > 150
            layer = 4;
        else
            layer = 2;
        end
    
    elseif age == 18 || age == 21 
        if cortical_layer_infor > 180
            layer = 4;
        else
            layer = 2;
        end   

    else  layer = NaN;  % Catch unhandled ages

    end
    
    % Convert image location string (e.g., 'L1') to numeric
    image_location =  regexprep(image_location_infor, '^L', '');
    
    %% Store extracted metadata
    all_name(end+1,1) = session_name;
    all_age (end+1,1) = age;
    all_LitterID(end+1,1) = Litter_ID;
    all_Layer(end+1,1) = layer;
    all_SubjectID(end+1,1) = subject_id;
    all_location(end+1,1) = image_location;
        
    %% Store clustering output results
    all_total_cell_number(end+1,1) = output_result.total_cell_number;
    all_NCl_beforeStat(end+1,1) = output_result.NCl_beforeStat;
    all_silhs_mean_beforeStat(end+1,1) = output_result.silhs_mean_beforeStat;
    all_NCl(end+1,1) = output_result.NCl;
    all_No_assemblies(end+1,1) =output_result.No_assemblies;
    all_S_assemblies(end+1,1) = output_result.S_assemblies;
    all_M_assemblies(end+1,1) = output_result.M_assemblies;
    all_Cells_not_in_assembly(end+1,1) = output_result.Cells_not_in_assembly;
    all_Cells_in_one_assembly(end+1,1) = output_result.Cells_in_one_assembly;
    all_Cells_in_many_assembly(end+1,1) = output_result.Cells_in_many_assembly;
    
    % Clear the variable for next loop
    clearvars output_result
end
%%
%% Define column names for output table
variable_names = {'LitterID', 'Layer','Age', 'SubjectID', 'Location', 'name',...
                  'total_cell_number', 'NCls_beforeStat','silhs_mean_beforeStat','NCls',...
                  'No_assemblies','S_assemblies','M_assemblies','Cells_not_in_assembly',...
                  'Cells_in_one_assembly','Cells_in_many_assembly'};

%% Combine all data into a summary table
tab = table(all_LitterID,all_Layer,all_age,all_SubjectID,all_location,all_name,...
            all_total_cell_number,all_NCl_beforeStat,all_silhs_mean_beforeStat,all_NCl,...
            all_No_assemblies,all_S_assemblies,all_M_assemblies,...
            all_Cells_not_in_assembly,all_Cells_in_one_assembly, all_Cells_in_many_assembly,...
            'VariableNames',variable_names); 
%%
%% Load subject sex information from external list and merge
subject_list_path = fullfile(root_path,"Subject_list.xlsx");
subject_list_table = readtable(subject_list_path,"VariableNamingRule","preserve");
subject_list_table = subject_list_table(:,1:2); % only keep "subjectID" and "sex" column

% Merge by SubjectID to append 'Sex' information
summary_tab = join(tab,subject_list_table,"LeftKeys","SubjectID","RightKeys","SubjectID");

%% Write final summary table to Excel           
write_file_path = fullfile(root_path,"Summary_file.xlsx");
writetable(summary_tab, write_file_path)
writetable(summary_tab, write_file_path, 'WriteVariableNames', true, 'WriteMode', 'overwritesheet')
disp("Finish summarizing data");