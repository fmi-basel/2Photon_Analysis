%% Collect Dataset for 2P Classification 
% (1) Use Masterstructure to assemble a Predictor Dataset, then (2) save the
% neuron information in a corresponding Matrix that allows to be loaded
% later for manual sorting and matched to the automatic Scoring
 
% written by Julian Hinz @ Luthi lab, 2021
if ~isempty(gcp('nocreate'))
    delete(gcp)
end

parpool(19);

clear all
close all
clc

% Start Directory
Direc = '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Imaging_2p';

% Add path that are important for Analysis
addpath(genpath('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Github\gluthi-2Photon'))

% Select Folders to Scan
Files = uipickfiles('FilterSpec', Direc);

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
        cell2mat(cellfun(@numel, {subFolders(:).name},'UniformOutput',false)) < 20);

    subFolders_all = subFolders(index_Ses);
    subFolders = sort_nat({subFolders_all.name});
    
    Path_to_Neural_dat = [Path_to_Neural_dat strcat([Filepath '\'], subFolders)];
end


% Now run through the neural data, either create the dataset for
% classification or load it, then save neural data in convenient Format

Matrix_Session = cell(size(Path_to_Neural_dat, 2), 1);
Neurons_Session = cell(size(Path_to_Neural_dat, 2), 1);
PointerAnimal_Session = cell(size(Path_to_Neural_dat, 2), 1);
Index_Count = 0;

for x = 1:size(Path_to_Neural_dat, 2)
    
    FilePath_cur = [Path_to_Neural_dat{x} '\processed_data_Folder'];
    
    if ~exist(FilePath_cur)
        disp('File has not been processed !')
        continue
    end
    
    if exist([FilePath_cur '\Results_CNMF.mat']) > 0
        % Extract Masterstructure
        load([FilePath_cur '\Timing_Information.mat'])
        load([FilePath_cur '\Corr_Imgs.mat'])
        load([FilePath_cur '\Results_CNMF.mat'])
        
        Master_structure = Extract_Components(neuron, obj, Timing_Info.FS);
    end

    % Check if an initial Manual Sorting has already taken place
    if exist([FilePath_cur '\Flagged_Components.mat']) > 0
        load([FilePath_cur '\Flagged_Components.mat'])
    else
        Refinement = [repelem({nan}, obj.PlanesImaged); repelem({nan}, ...
            obj.PlanesImaged)]';
    end
    
    % Now assemble the 2 required sets: Lets start with the predictor set
    Matrix_Cellspace_tmp = cell(obj.PlanesImaged, 1);
    Cell_Info_tmp = cell(obj.PlanesImaged, 1);
    PointerAnimal = cell(obj.PlanesImaged, 1);
    Neurons_results = cell(obj.PlanesImaged, 1);
    
    for xx = 1:obj.PlanesImaged
        
        % Assemble predictor structure
        Cur_Plane = Master_structure.(['Plane_' num2str(xx)]);
        
        Size = size(Cur_Plane.Corr_imgs);
        Distance_from_center = pdist2([Size(1)/2 Size(2)/2], ...
            Cur_Plane.Centroids) / Size(1);
        Index = [Index_Count+1:Index_Count+ size(Distance_from_center, 2)]';
        
        if all(~isnan(Refinement{xx, 1}))
            Decision = Refinement{xx, 1};
        else
            Decision = nan(size(Index, 1), 1);
        end
        
        Matrix_Cellspace_tmp{xx} = [Index Cur_Plane.Mean_2_std Cur_Plane.Extracted_values' ...
            Cur_Plane.Area Cur_Plane.circularities Cur_Plane.Perimeter ...
            Cur_Plane.Eccentricity Distance_from_center' ...
            Cur_Plane.Std2crossings Cur_Plane.Correlations' Cur_Plane.AbsoluteSpikeCount ...
            Cur_Plane.Dec_Rawcorr' Cur_Plane.CorrImgs2MedianCorr' Cur_Plane.Sizerel2max ...
            Cur_Plane.Corr_2Local Cur_Plane.Overlaps Decision];
        
        % Now assemble the Data Structure containing Informations regarding 
        % all neurons matched to indices
        
        CurResults = neuron.(['Plane_' num2str(xx)]);
        
        % Lets start with spatial Footprint
        if size(CurResults.Cn, 1) > 256
            Neurons_Spatial = cellfun(@(d) imresize(d, [256, 256], 'bicubic'),...
                cellfun(@(c) reshape(c, size(CurResults.Cn, 1), size(CurResults.Cn, 2)), ...
                mat2cell(CurResults.A', ones(size(CurResults.A, 2), 1), ...
                size(CurResults.A, 1)), 'UniformOutput', false), 'UniformOutput', false);
            
            Neurons_Cen = cellfun(@(c) round(c/(size(CurResults.Cn, 1)/256)), ...
                mat2cell(Cur_Plane.Centroids, ones(size(Cur_Plane.Centroids, 1), 1), ...
                size(Cur_Plane.Centroids, 2)), 'UniformOutput', false);  
        else
            Neurons_Spatial = cellfun(@(c) reshape(c, size(CurResults.Cn, 1), size(CurResults.Cn, 2)), ...
                mat2cell(CurResults.A', ones(size(CurResults.A, 2), 1), ...
                size(CurResults.A, 1)), 'UniformOutput', false);
            Neurons_Cen = mat2cell(Cur_Plane.Centroids, ones(size(Cur_Plane.Centroids, 1), 1), ...
                size(Cur_Plane.Centroids, 2));  
        end
        
        Neurons_Deconv = mat2cell(Cur_Plane.C_dec, ones(size(Cur_Plane.C_dec, 1), 1), ...
                size(Cur_Plane.C_dec, 2));     
        Neurons_Fdff = mat2cell(Cur_Plane.F_dff, ones(size(Cur_Plane.F_dff, 1), 1), ...
                size(Cur_Plane.F_dff, 2));    
        Neurons_Raw = mat2cell(CurResults.C, ones(size(CurResults.C, 1), 1), ...
            size(CurResults.C, 2));  
        
        Neurons_Timing = repelem({Timing_Info.FS}, size(Neurons_Raw, 1))';
%         Neurons_CN = repelem({obj.(['Cor_Images_' num2str(xx)]).Corr_Img(:,:,4)}, ...
%             size(Neurons_Raw, 1))';
         
        Neurons_results{xx} = [Neurons_Spatial Neurons_Deconv Neurons_Fdff ...
            Neurons_Cen Neurons_Raw Neurons_Timing]; % Neurons_CN
        
        % Pointer to Plane and Folder
        PointerAnimal{xx} = [mat2cell(Index, ones(size(Neurons_Raw, 1), 1), 1) ...
            repelem({FilePath_cur}, size(Neurons_Raw, 1))' ...
            repelem({['Plane_' num2str(xx)]}, size(Neurons_Raw, 1))'];
        
        % Increase Index to track the predictor to neuron relationship
        Index_Count = Index_Count + size(Distance_from_center, 2);
        
    end
    
    % Collection of Session Information
    Matrix_Session{x} = cat(1, Matrix_Cellspace_tmp{:});
    
    % Neuron Information
    Neurons_Session{x} = cat(1, Neurons_results{:});
    
    % Pointer to animal folder and Plane
    PointerAnimal_Session{x} = cat(1, PointerAnimal{:});
    
    disp(['Finished with Folder: ' Path_to_Neural_dat{x}])
    
end

% Now Save the Data
CompleteMatrix = cat(1, Matrix_Session{:});
savefast('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Classification\TwoPhoton\ClassMatrix.mat', 'CompleteMatrix');

CompleteNeurons = cat(1, Neurons_Session{:});
savefast('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Classification\TwoPhoton\Neurons.mat', 'CompleteNeurons');

CompletePointer = cat(1, PointerAnimal_Session{:});
savefast('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Classification\TwoPhoton\Pointers.mat', 'CompletePointer');







