%% Generate Heuristics PostHoc

clear
close
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the extraction
Run_Extraction = true;

% Start Directory
Direc = '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Imaging_2p';

% Add path that are important for Analysis
addpath(genpath('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Github\gluthi-2Photon'))


Files = uipickfiles('FilterSpec', Direc);
One_up = strfind(Files, '\');

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Run_Extraction
    for x = 1:size(Path_to_Neural_dat, 2)
        try
            load([Path_to_Neural_dat{x} '\processed_data_Folder\Results_CNMF.mat'])
            load([Path_to_Neural_dat{x} '\processed_data_Folder\Corr_Imgs.mat'])
            load([Path_to_Neural_dat{x} '\processed_data_Folder\Timing_Information.mat'])
            
            Master_structure = Extract_Components(neuron, obj, Timing_Info.FS);
            savefast([Path_to_Neural_dat{x} '\processed_data_Folder\Master_structure.mat'], 'Master_structure')
        catch
            disp('Couldnt load Folder')
            continue
        end
        
        disp(['Finished Folder ' num2str(x) ' out of ' num2str(size(Path_to_Neural_dat, 2)) ' Folders !'])
    end
end

GoodTraces = cell(size(Path_to_Neural_dat, 2), 1);
BadTraces = cell(size(Path_to_Neural_dat, 2), 1);
Predictors = cell(size(Path_to_Neural_dat, 2), 1);

for x = 1:size(Path_to_Neural_dat, 2)
    
    % Load data
    try
        load([Path_to_Neural_dat{x} '\processed_data_Folder\Master_structure.mat'])
        load([Path_to_Neural_dat{x} '\processed_data_Folder\Flagged_Components.mat'])
    catch
        disp('Couldnt load Folder')
        continue
    end
    
    % Extract individual Planes
    Planes2Runthrough = fieldnames(Master_structure);
    GoodTraces_tmp = cell(size(Planes2Runthrough, 1), 1);
    BadTraces_tmp = cell(size(Planes2Runthrough, 1), 1);
    Predictors_tmp = cell(size(Planes2Runthrough, 1), 1);
    
    for xx = 1:size(Planes2Runthrough, 1) 
        
        % Extract Traces
        GoodTraces_tmp{xx} = Master_structure.(Planes2Runthrough{xx}).F_dff(~Refinement{xx, 2}, :);
        BadTraces_tmp{xx} = Master_structure.(Planes2Runthrough{xx}).F_dff(Refinement{xx, 2}, :);
        
        % Now extract Predictors
        tmp = Master_structure.(Planes2Runthrough{xx});
        Predictors_tmp{xx} = [tmp.Correlations' tmp.Eccentricity tmp.Perimeter ...
            tmp.circularities tmp.Area tmp.Mean_2_std tmp.Extracted_values' ...
            tmp.distance_from_centre' tmp.AbsoluteSpikeCount tmp.Dec_Rawcorr' ...
            tmp.CorrImgs2MedianCorr' tmp.Sizerel2max tmp.Corr_2Local tmp.Overlaps ...
            tmp.Std2crossings Refinement{xx, 2}];
    end
    
    GoodTraces{x} = cat(1, GoodTraces_tmp{:});
    BadTraces{x} = cat(1, BadTraces_tmp{:});
    Predictors{x} = cat(1, Predictors_tmp{:});
end

% Subsample predictors
rng(1);
Predictors_conc = cat(1, Predictors{:});
Indices = datasample(1:size(Predictors_conc, 1), 15000, 'Replace',false);
Predictors_Tianlin = Predictors_conc(Indices, :);
Predictors_conc(Indices, :) = [];

% Save results
savefast('Z:\hinzjuli\Results\TianlinTest\Predictors.mat', 'Predictors_Tianlin')
savefast('Z:\hinzjuli\Results\TianlinTest\Predictors_Test.mat', 'Predictors_conc')
savefast('Z:\hinzjuli\Results\TianlinTest\BadTraces.mat', 'BadTraces')
savefast('Z:\hinzjuli\Results\TianlinTest\GoodTraces.mat', 'GoodTraces')
