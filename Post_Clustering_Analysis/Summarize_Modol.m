function Summarize_Modol(root_path)
% Summarize outputs from Modol code
% 
% First, let's turn off a warning about loading parallel clusters; this happens because the Modol code saves all variables including the information about parallel workers, but this information cannot be loaded back in.
warning('off', 'parallel:cluster:CannotLoadCorrectly')
% 
% Set the path of the litter
% root_path = "G:\Michael\20210129-0208_jRGECO1a_Chrm3q-Ncre19-2";

% Create containers for the variables of interest
ages = ["P11", "P13", "P15", "P18", "P21", "P30"];

NCls = [];
silhs = [];
S_assemblies = [];
M_assemblies = [];
names = [];

% Loop through all the animal folders in the root path, load the modol outputs generated for that animal, and then add the average sihlouette value, the number of significant cell assemblies, and the number of S/M-assembly GCEs
age_dir = dir(root_path);
%disp("Age dir:")
for x = 1:length(age_dir)
    disp(age_dir(x).name)
    % Only loop through age folders
    if strcmp(age_dir(x).name, ".") || strcmp(age_dir(x).name, "..") || isempty(regexp(age_dir(x).name, 'P[0-9]{2}', 'once'))
        continue
    end
    % Get the age of interest
    age_idx = regexp(age_dir(x).name, 'P[0-9]{2}');
    age_name = age_dir(x).name;
    age_name = age_name(age_idx:age_idx+2);
    ages_idx = find(contains(ages, age_name));
    age = ages(ages_idx);
    
    % Loop through age directory
    anim_dir = dir(root_path+"/"+age_dir(x).name);
    disp(anim_dir)
    for jdx = 1:length(anim_dir)
        disp(anim_dir(jdx).name)
        % Skip stimulation data
        if ~isempty(strfind(anim_dir(jdx).name, "stim")) || strcmp(anim_dir(jdx).name, ".") || strcmp(anim_dir(jdx).name, "..")
            continue
        elseif isempty(regexp(anim_dir(jdx).name, 'P[0-9]{2}', 'once')) && isempty(regexp(anim_dir(jdx).name, '#[0-9]', 'once'))
            continue 
        end
        
        % Load Modol outputs
        load_loc = root_path+"/"+age_dir(x).name+"/"+anim_dir(jdx).name+"/Modol_outputs/all.mat";
        load_loc = fullfile(load_loc);
        disp("Load loc:")
        disp(load_loc)
        try
            load(load_loc);
        catch ME
            fprintf('%s could not be loaded', load_loc)
            continue 
        end
        % If age not in name, add it
        if isempty(regexp(anim_dir(jdx).name, 'P[0-9]{2}', 'once'))
            name = strcat(age, "_", anim_dir(jdx).name);
        else
            name = anim_dir(jdx).name;
        end
        % Total number of cells
        total = size(detected_events, 1) - length(A0);
        % Stack variables into containers
        names = [names; string(name)];
        NCls = [NCls; NCl];
        silhs = [silhs; mean(sCl)];
        S_assemblies = [S_assemblies; length(A1)/total*100];
        M_assemblies = [M_assemblies; length(A2)/total*100];
    end
end


% Save the summarized outputs (in the root directory/Summary_Modol)
if ~exist(root_path+"Summary_Modol/", 'dir')
    mkdir(root_path+"Summary_Modol/")
end
tab = rows2vars(table(names, NCls, silhs, S_assemblies, M_assemblies));
writetable(tab, root_path+"Summary_Modol/Summary_Modol.xlsx", 'WriteVariableNames', 0, 'WriteMode', 'overwritesheet')
