%% 2P Imaging Analysis
% Julian Hinz @ Luthi lab 2018 - 2025
% To use this script, please refer to the pdf detailing the
% required Folder directory structure and naming conventions

if ~isempty(gcp('nocreate'))
    delete(gcp)
end 

clear
close
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Direc = '';

% Specification for the Image analysis
min_Size_Components = 10;        
T_DS_Factor_standard = 1; % Temporal downsampling factor - Check lower to change
S_DS_Factor = 1; % spatial downsampling factor
vis = false; % visualization during analysis
FS = 10;
DS_Factor_vis = FS*2;

% Do you want memory mapping
MemMap = true;

% Scanning mode : Bidirectional or SingleDirection
Scanning_mode = 'Bidirectional'; % Bidirectional

% Limit number of Cores
last_comp_threads = maxNumCompThreads(4);

% Max Runtime MC and CNMF
Max_Runtime = 5;
parpoolTest = 2;

% Number of parallel workers
num_workers_t = 5;

% Square processing
square = true;

% Classification through SVM
SVM_classification = true;

% Post Processing
RegionSub = true;
manual_Selection = true;

% Add folder directory automatically
filePath = matlab.desktop.editor.getActiveFilename;
One_up = strfind(filePath, '\');
addpath(genpath(filePath(1:One_up(end))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select files and generate experimental details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Files = uipickfiles('FilterSpec', Direc);
One_up = strfind(Files, '\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Motion Correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify paths to the neural data
Path_to_Neural_dat = [];

for i = 1:size(Files, 2)

    % Current Folder    
    Filepath = Files{i};

    % Automatically detect folders that need to be processed
    files = dir(Filepath);% Get a list of all files and folders in this folder.
    dirFlags = [files.isdir]; % Get a logical vector that tells which is a directory.
    subFolders = files(dirFlags); % Extract only those that are directories.
    
    % Sort the folders and select the correct ones based on length and name
    index_Ses = find(cellfun(@(x) any(strfind(x, 'Session')), {subFolders(:).name}) & ...
        cell2mat(cellfun(@numel, {subFolders(:).name}, 'UniformOutput', false)) < 20);
    subFolders_all = subFolders(index_Ses);          
    subFolders = sort_nat({subFolders_all.name});

    Path_to_Neural_dat = [Path_to_Neural_dat strcat([Filepath '\'], subFolders)];

end

% Now loop over data and each time offer selection to reprocess part of the
% data - Currently no possible to select parts of the data

Run_time = 0;
excluded_all = cell(3, 1);

while Run_time <= 3 && ~isempty(Path_to_Neural_dat)
    
    PathExcluded = zeros(size(Path_to_Neural_dat, 2), 1);
    Run_time = Run_time + 1;
    
    for i = 1:size(Path_to_Neural_dat, 2)
        
        % Reset number of workers
        num_workers = num_workers_t;
        
        % Reset Downsampling
        T_DS_Factor = T_DS_Factor_standard;
        
        % Current Folder    
        Filepath = Path_to_Neural_dat{i};
        if ~exist([Filepath '\processed_data_Folder'])
            mkdir([Filepath '\processed_data_Folder'])
        elseif exist([Filepath '\processed_data_Folder\Corr_Imgs_' num2str(Run_time) '.mat']) ~= 0 || ...
            exist([Filepath '\processed_data_Folder\Corr_Imgs.mat']) ~= 0

            fprintf('\n\n###################################################################\n')
            fprintf(' Session was already processed, skip !!!\n')
            fprintf('###################################################################\n\n')
            continue
        end

        Loading = false;
        counter = 1;
        while ~Loading               
            try           
                % Concatenation of tif Files
                % Insert an automatic mechanism to archive - Then save the path
                Patch_size = 10000;
                [num_planes, Img_Folder, Length_Recording, no_imgs] = Stack_Tiff_Bruker([Filepath '\'], ...
                    Patch_size, MemMap, Scanning_mode, num_workers, square);
                Loading = true;
            catch
                counter = counter + 1;
                Directories = dir(fullfile([Filepath '\processed_data_Folder\'], '*.mat'));   
                cellfun(@delete, (strcat([Filepath '\processed_data_Folder\'], ...
                {Directories(cellfun(@(x) any(strfind(x, 'Raw')), {Directories.name})).name})));
                fprintf('\n\n###################################################################\n')
                fprintf(' Error Loading - Retrying !!\n')
                fprintf('###################################################################\n\n')
            end
            
            if counter >= 2 && ~Loading  
                fprintf('\n\n###################################################################\n')
                fprintf(' Waiting 60 min before retrying - Excluding Network issues ! \n')
                fprintf('###################################################################\n\n')
                pause(3600)
            end
            
            if counter > Max_Runtime && ~Loading  
                fprintf('\n\n###################################################################\n')
                fprintf(' Maximum Runtime exceeded - Exiting analysis !!\n')
                fprintf('###################################################################\n\n')
                fprintf('\n')
                
                
                error(['Loading failed to process Files in Folder: ' Filepath ...
                     ' - Please check whats up with this file !']) 
                fprintf('\n')
            end
        end
        
        if no_imgs
            PathExcluded(i) = 1;
            continue
        end
        
        % Imaging Analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Motion Correction

        % Normcorre
        % see https://github.com/flatironinstitute/NoRMCorre
        % see https://www.sciencedirect.com/science/article/pii/S0165027017302753?via%3Dihub
        % Rigid Motion correction of the 2P data, non-Rigid MC is not yet
        % implemented, but is available in the same toolbox
        
        if Run_time == 1
            Planes_to_process = 1:num_planes;
            Selection = [];
        else
            Planes_to_process = find(quality_passed{i} == 0);
            Selection = Selection_all{i};
        end
                
        % There were some problems with the MC, not being able to write,
        % which could be due to glitches in the network. Therefore its
        % wrapped in a try catch statement, which will repeat it a couple
        % of times - Will error out, if its too often then IT needs to see
        % how to fix it
        
        MC_completed = false;
        counter = 1;
        while ~MC_completed               
            try
                [Movies, Templates_rg, Shifts_summarized] = TwoP_Image_Registration(Filepath, ...
                    MemMap, T_DS_Factor, S_DS_Factor, DS_Factor_vis, Planes_to_process, ...
                    Selection, Patch_size, Run_time, num_workers);      
                MC_completed = true;
            catch
                counter = counter + 1;
                fprintf('\n\n###################################################################\n')
                fprintf(' Error Motion Correcting - Retrying !!\n')
                fprintf('###################################################################\n\n')
            end
            
            if counter >= 2
                fprintf('\n\n###################################################################\n')
                fprintf(' Waiting 60 min before retrying - Excluding Network issues ! \n')
                fprintf('###################################################################\n\n')
                num_workers = parpoolTest;
                pause(3600)
            end
            
            if counter > Max_Runtime 
                fprintf('\n\n###################################################################\n')
                fprintf(' Maximum Runtime exceeded - Exiting analysis !!\n')
                fprintf('###################################################################\n\n')
                fprintf('\n')
                                
                error(['NormCorre failed to process Files in Folder: ' Filepath ...
                     ' - Please check whats up with this file !']) 
                fprintf('\n')
            end
        end
    
        if ~isempty(Movies)
            savefast([Filepath '\processed_data_Folder\Movies_' num2str(Run_time) '.mat'] , 'Movies');
            savefast([Filepath '\processed_data_Folder\Motion_correction_templates_' num2str(Run_time) '.mat'] , 'Templates_rg');
            savefast([Filepath '\processed_data_Folder\Shifts_' num2str(Run_time) '.mat'] , 'Shifts_summarized');
        end
        
        clear Movies Templates_rg Shifts_summarized
        
        %%% Calculate the Correlation Image
        % This is what I use to exctract components later
        
        [obj] = Corr_Img_Segmentation(Filepath, vis, num_planes);  
        obj.PlanesImaged = num_planes;
        obj.T_DS_Factor = T_DS_Factor;
        obj.S_DS_Factor = S_DS_Factor;
        obj.DS_Factor_vis = DS_Factor_vis;
        
        savefast([Filepath '\processed_data_Folder\Corr_Imgs_' num2str(Run_time) '.mat'], 'obj')
        clear obj  
        
        fprintf('\n\n###################################################################\n')
        disp([' Finished Motion Correction of Folder: ' Filepath])
        fprintf('###################################################################\n\n')
    end
    
    Path_to_Neural_dat = Path_to_Neural_dat(~logical(PathExcluded));
    
    % Load all movies, in order not to clock memory during processing
    Movies2Check_tmp = cell(size(Path_to_Neural_dat, 2), 1);
    objs_coll = cell(size(Path_to_Neural_dat, 2), 1);
    delete_var = false(size(Path_to_Neural_dat, 2), 1);
    Excluded = cell(size(Path_to_Neural_dat, 2), 1);
    
    for i = 1:size(Path_to_Neural_dat, 2)
        if exist([Path_to_Neural_dat{i} '\processed_data_Folder\Movies.mat']) ~= 0
            
            Qualitys = dir([Path_to_Neural_dat{i} '\processed_data_Folder\QualityCheck*']);
            delete_var(i) = true;
            if size(Qualitys, 1) > 1
                [~, idx] = sort_nat({Qualitys.name}, 'descend');
                load([Path_to_Neural_dat{i} '\processed_data_Folder\' Qualitys(idx(1)).name]);
                if ~all(quality_passed_Session)
                    Excluded{i} = Path_to_Neural_dat{i};
                end
            end  
            continue
        else
            Movies2Check_tmp{i} = load([Path_to_Neural_dat{i} '\processed_data_Folder\Movies_' num2str(Run_time) '.mat']);
            objs_coll{i} = load([Path_to_Neural_dat{i} '\processed_data_Folder\Corr_Imgs_' num2str(Run_time) '.mat']);
        end
    end
            
    Path_to_Neural_dat(delete_var) = [];
    objs_coll(delete_var) = [];
    Movies2Check_tmp(delete_var) = [];
    
    % Quality Control for the Motion Correction
    MC_Selection = cell(size(Path_to_Neural_dat, 2), 1);
    quality_passed = cell(size(Path_to_Neural_dat, 2), 1);
    Selection_all = cell(size(Path_to_Neural_dat, 2), 1);
    
    for i = 1:size(Path_to_Neural_dat, 2)

        % Current Folder    
        Filepath = Path_to_Neural_dat{i};
        obj = objs_coll{i}.obj;
        quality_passed_Session = zeros(obj.PlanesImaged, 1);
        MC_Selection_tmp = cell(obj.PlanesImaged, 1);
                
        fig = figure('position', [1300 100 500 500]); 
        for z = 1:obj.PlanesImaged
            if obj.PlanesImaged > 1
                subplot(ceil(obj.PlanesImaged/3), 3, z)
                imagesc(obj.(['Cor_Images_' num2str(z)]).Corr_Img(:,:,4)); axis off; colormap bone; box off; 
            else
                imagesc(obj.(['Cor_Images_' num2str(z)]).Corr_Img(:,:,4)); axis off; colormap bone; box off; 
            end
        end
        
        One_up = strfind(Filepath, '\');  
        if obj.PlanesImaged > 1
            suptitle(Filepath(One_up(end)+1:end))
        else
           title(Filepath(One_up(end)+1:end)) 
        end
        
        Handle = implay(Movies2Check_tmp{i}.Movies, 60);
        Handle.Parent.Position = [100, 50, 1300, 1300];

        fprintf('Video okay ? Select the ones that are flawed by number, or d to delete session:    ');
        temp = input('', 's');
        
        if isempty(temp) 
            quality_passed_Session(:) = 1; 
            Selection_all{i} = [];
        elseif strcmpi(temp, 'd') 
            quality_passed_Session(:) = 1; 
            Excluded{i} = Filepath;
            Selection_all{i} = [];
        else
            try  
                Numbers = cellfun(@str2num, regexp(temp, ['[1-' ...
                    obj.PlanesImaged ']'], 'match'));
            catch
                disp('No numbers selected, all numbers are selected !')
                Numbers = 1:obj.PlanesImaged;
            end
            
            [Selection_all{i}] = Define_ROIS_2P(obj.(['Cor_Images_' num2str(1)]).Corr_Img(:,:,4), 1);
            quality_passed_Session(:) = 1;
            quality_passed_Session(Numbers) = 0;            
        end
        
        try
            close(Handle)
        catch
            disp('Please dont close figures !')
        end
        
        try
            close(fig)  
        catch
            disp('Please dont close figures !')
        end
        
        savefast([Filepath '\processed_data_Folder\QualityCheck_' ...
        num2str(Run_time) '.mat'], 'quality_passed_Session');
        
        % Take care of the individual Raw files and delete them if they
        % were all properly processed
        
        if all(quality_passed_Session)
            Directories = dir(fullfile([Filepath '\processed_data_Folder\'], '*.mat'));            
            cellfun(@delete, (strcat([Filepath '\processed_data_Folder\'], ...
                {Directories(cellfun(@(x) any(strfind(x, 'Raw_')), {Directories.name})).name})));
            
            % Delete the numbering and rename all the variables, so that
            % there is no confusion with the numbering
            movefile([Filepath '\processed_data_Folder\Shifts_' num2str(Run_time) '.mat'], ...
                [Filepath '\processed_data_Folder\Shifts.mat'], 'f');
            movefile([Filepath '\processed_data_Folder\Corr_Imgs_' num2str(Run_time) '.mat'], ...
                [Filepath '\processed_data_Folder\Corr_Imgs.mat'], 'f');            
            movefile([Filepath '\processed_data_Folder\Movies_' num2str(Run_time) '.mat'], ...
                [Filepath '\processed_data_Folder\Movies.mat'], 'f');           
            movefile([Filepath '\processed_data_Folder\QualityCheck_' num2str(Run_time) '.mat'], ...
                [Filepath '\processed_data_Folder\QualityCheck.mat'], 'f');
            movefile([Filepath '\processed_data_Folder\Motion_correction_templates_' num2str(Run_time) '.mat'], ...
                [Filepath '\processed_data_Folder\Motion_correction_templates.mat'], 'f');    

            if Run_time > 1
                for b = 1:Run_time - 1
                    delete([Filepath '\processed_data_Folder\Shifts_' num2str(b) '.mat']); 
                    delete([Filepath '\processed_data_Folder\Corr_Imgs_' num2str(b) '.mat']);
                    delete([Filepath '\processed_data_Folder\Movies_' num2str(b) '.mat']); 
                    delete([Filepath '\processed_data_Folder\Motion_correction_templates_' num2str(b) '.mat']); 
                end
            end
        end
                
        % If quality of all is passed, delete all individual raw files
        MC_Selection_tmp{i} = MC_Selection_tmp;
        quality_passed{i} = quality_passed_Session;
    end
    
    clear Movies2Check_tmp
    fprintf('\n\n###################################################################\n')
    disp(['Finished Run ' num2str(Run_time) ', in total ' num2str(sum(cellfun(@all, quality_passed))) ...
        ' out of ' num2str(numel(Path_to_Neural_dat)) ' Sessions were successfully processed.'] )
    fprintf('###################################################################\n\n')
    
    excluded_all{Run_time} = Excluded(~cellfun(@isempty, Excluded));
    Path_to_Neural_dat(cellfun(@all, quality_passed)) = [];
    quality_passed(cellfun(@all, quality_passed)) = [];
    Selection_all(cellfun(@isempty, Selection_all)) = [];    
end

fprintf('\n\n###################################################################\n')
disp('Finished Motion Correction of all Folders, now entering Extraction !')
fprintf('###################################################################\n\n')

% Combine the excluded all file, so that excluded sessions can be skipped
% later

excluded_all = excluded_all(~cellfun(@isempty, excluded_all));
excluded_all = cat(1, excluded_all{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Segment Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Automatic Segmentation of Images with the CNMF
% Need to test whether the Session was cancelled somehow with the filepath
% saved in cell excluded

% Start parallel pool
num_workers = num_workers_t;
parpool(num_workers);

for i = 1:size(Files, 2)
    
    % Current Folder    
    Filepath = Files{i};
    
    % Automatically detect folders that need to be processed
    files = dir(Filepath);% Get a list of all files and folders in this folder.
    dirFlags = [files.isdir]; % Get a logical vector that tells which is a directory.
    subFolders = files(dirFlags); % Extract only those that are directories.
    

    % Sort the folders and select the correct ones based on length and name
    index_Ses = find(cellfun(@(x) any(strfind(x, 'Session')), {subFolders(:).name}) & ...
        cell2mat(cellfun(@numel, {subFolders(:).name},'UniformOutput',false)) < 20);

    subFolders_all = subFolders(index_Ses);
    subFolders = sort_nat({subFolders_all.name});
    
    for ii = 1:size(subFolders, 2)
                
        if any(ismember([Filepath '\' subFolders{ii} '\processed_data_Folder'], excluded_all))
            disp('File was excluded from further analysis !')
            continue
        end
        
        if num_workers == parpoolTest
            if ~isempty(gcp('nocreate'))
                delete(gcp)
            end
            parpool(num_workers_t);
        end
        
        Extraction = false;
        counter = 1;
        while ~Extraction % Run extraction three times
            try
                Perform_CNMF2([Filepath '\' subFolders{ii} '\processed_data_Folder'], Max_Runtime, FS)
                Extraction = true;
            catch
                disp('File not processed, trying again !')
                if counter == 1 % first decrease number of cores as likely cause of crash
                    if ~isempty(gcp('nocreate'))
                        delete(gcp)
                    end
                    parpool(parpoolTest);
                    % if the number of cores and also network problems are 
                    % excluded through rerunning a copuple of times, exit 
                    % and display error message
                elseif counter > Max_Runtime 
                   disp(['File not processed, please check: ' Filepath '\' ...
                       subFolders{ii} '\processed_data_Folder'])
                   Extraction = true;
                end
                counter = counter + 1;
            end
        end
        
        fprintf('\n\n###################################################################\n')
        disp([' Finished Extraction of Folder: ' Filepath '\' subFolders{ii}])
        fprintf('###################################################################\n\n')       
    end
end

delete(gcp)

fprintf('\n\n###################################################################\n')
disp('Finished Extraction of all Folders, now entering Post-Processing !')
fprintf('###################################################################\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Curation of Segmented Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, automatic classification, then manual curation of the rest of the
% data - Using SVM classifier

ReProcess_flag = false;

if SVM_classification
    
    for i = 1:size(Files, 2)

        % Current Folder    
        Filepath = Files{i};

        % Automatically detect folders that need to be processed
        files = dir(Filepath);% Get a list of all files and folders in this folder.
        dirFlags = [files.isdir]; % Get a logical vector that tells which is a directory.
        subFolders = files(dirFlags); % Extract only those that are directories.


        % Sort the folders and select the correct ones based on length and name
        index_Ses = find(cellfun(@(x) any(strfind(x, 'Session')), {subFolders(:).name}) & ...
            cell2mat(cellfun(@numel, {subFolders(:).name},'UniformOutput',false)) < 20);

        subFolders_all = subFolders(index_Ses);
        subFolders = sort_nat({subFolders_all.name});

    end
    for ii = 1:size(subFolders, 2)
        if ~exist([Filepath '\' subFolders{ii} '\processed_data_Folder\'...
                'PospP_CNMF_SVM.mat']) || ReProcess_flag
            Automatic_Classification([Filepath '\' subFolders{ii} '\processed_data_Folder'])
        end

        fprintf('\n\n###################################################################\n')
        disp([' Finished Classification of Folder: ' Filepath '\' subFolders{ii}])
        fprintf('###################################################################\n\n')
    end
else

    for i = 1:size(Files, 2)

        % Current Folder    
        Filepath = Files{i};

        % Automatically detect folders that need to be processed
        files = dir(Filepath);% Get a list of all files and folders in this folder.
        dirFlags = [files.isdir]; % Get a logical vector that tells which is a directory.
        subFolders = files(dirFlags); % Extract only those that are directories.

        % Sort the folders and select the correct ones based on length and name
        index_Ses = find(cellfun(@(x) any(strfind(x, 'Session')), {subFolders(:).name}) & ...
            cell2mat(cellfun(@numel, {subFolders(:).name},'UniformOutput',false)) < 20);

        subFolders_all = subFolders(index_Ses);
        subFolders = sort_nat({subFolders_all.name});

  
    end
    for ii = 1:size(subFolders, 2)
        CNMF_refinement([Filepath '\' subFolders{ii} '\processed_data_Folder'], ...
            RegionSub, manual_Selection) 

        fprintf('\n\n###################################################################\n')
        disp([' Finished Extraction of Folder: ' Filepath '\' subFolders{ii}])
        fprintf('###################################################################\n\n')

    end
      
end




