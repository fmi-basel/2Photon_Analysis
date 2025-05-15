%% Seperation of files and subsequent zipping -> Needs to be implemented
% Author: Julian Hinz, 2018 (julian.hinz@fmi.ch)
% Function to seperate Individual Tiff files into several seperated tif
% stacks

function [num_planes, Files_Img, Length_Recording, no_imgs] = Stack_Tiff_Bruker(...
    Files, Patch_size, MemMap, Scanning_mode, num_workers, ...
    square)

    % The function expects following inputs:
    % Files = path to parent folder with individual tiff files
    % num_channels = number of channels that were used (Integer)
    % Patch_size = to be memory efficent the Tif-files are only loaded sequentially
    % usually sizes of 5000 should be fine
    % Statement if memory mapping is used or not, default true

    % Set options for saving
    options.compress = 'no';
    options.message = false;
    options.append = true;
    options.overwrite = false;
    options.big = true;
        
    if MemMap
        filet = '*.mat';        
    else
        filet = '*.tif';
    end

    Files = [Files 'Imaging\'];

    % Find the .tif files in specified Folder
    Temporary = dir(Files);
    Temporary = Temporary(cellfun(@numel, {Temporary.name}) > 10);

    if size(Temporary, 1) == 1
        Files_Img = [Files Temporary.name '\'];
    elseif size(Temporary, 1) > 1 && size(Temporary, 1) < 1000 
        error('Check whats weird here !')  
    else
        Files_Img = Files; 
    end

    % Find pos of processed data folder
    One_up = strfind(Files, '\');
    searchString2 = fullfile([Files(1:One_up(end-1)) '\processed_data_Folder\'], filet);
    loc = One_up(end-1);
   
    searchString = fullfile(Files_Img, '*.tif');
    Filenames = dir(searchString);

    if isempty(Filenames)
        no_imgs = true;
        num_planes = nan;
        Length_Recording = nan;
        disp('Terminating function no img files !')
        return
    else
        no_imgs = false;
    end
    
    % Read how many planes were recorded from filename  
    One_up_Fold = strfind(Files_Img, '\');
    FolderName = Files_Img(One_up_Fold(end-1)+1:end-1);
    Pos = strfind(lower(FolderName) , 'planes'); 
    Pos2 = strfind(FolderName , '_'); 
    Indices = Pos2 - Pos;
    PotNum = str2double(FolderName(Pos + min(Indices(Indices > 0)) - 1));   

    if isempty(PotNum) || isnan(PotNum)
        num_planes = str2double(Filenames(1).name(Pos + Indices(find(Indices == ...
             min(Indices(Indices > 0))) + 1) - 1));
    elseif isnumeric(PotNum)
        num_planes = PotNum;
    else
        error('Couldnt read number of Planes, please adjust script !')  
    end
    
    % Check if multiplane aquisition was terminated prematurely
    if size(Filenames, 1) > 50
        if mod(size(Filenames, 1), num_planes) ~= 0
            Subtract = mod(size(Filenames, 1), num_planes);
        else
            Subtract = 0;
        end
 
        Length_Recording = (size(Filenames, 1) - Subtract)/num_planes;
        NewStacks = false;
    else
        InfoImgs = nan(size(Filenames, 1), 3);
        for k = 1:size(Filenames, 1)
             InfoTif = imfinfo([Filenames(k).folder '\' Filenames(k).name]);
             InfoImgs(k, :) = [InfoTif(1).Width InfoTif(1).Height size(InfoTif, 1)];
        end
        Length_Recording = sum(InfoImgs(:, 3));
        NewStacks = true;
    end
    
    % Find folder up one level 
    Filenames2 = dir(searchString2);
    index2 = cellfun(@(x) any(strfind(x, 'Raw')), {Filenames2.name});
    
    if ~isempty(index2)
        if any(index2 ~= 0)
            disp('Files have already been stacked, terminating the function')
            return
        end
    end

    if mod(num_planes, Patch_size) ~= 0 && num_planes > 1
        Patch_size = Patch_size - mod(Patch_size, num_planes);
    end

    Iterations = ceil(Length_Recording/Patch_size);
    
    if MemMap        
        im_info = imfinfo([Filenames(1).folder '\' Filenames(1).name]);
        
        if im_info(1).Height ~= im_info(1).Width && square
            MinimumAxis = min([im_info(1).Height im_info(1).Width]);
            X_Axis = MinimumAxis;
            Y_Axis = MinimumAxis;
            Resizing = [im_info(1).Height/MinimumAxis im_info(1).Width/MinimumAxis];
        else
            X_Axis = im_info(1).Height;
            Y_Axis = im_info(1).Width;
            Resizing = [1 1];
        end
        
        if num_planes > 1
            Y = zeros(X_Axis, Y_Axis, Length_Recording, num_planes, 'uint16');
        else
            Y = zeros(X_Axis, Y_Axis, Length_Recording, 'uint16');
        end
        
        savefast([Files(1:loc) '\processed_data_Folder\Raw.mat'], 'Y');
    
        Path1 = [Files(1:loc) '\processed_data_Folder\MC_r.mat'];        
        copyfile([Files(1:loc) '\processed_data_Folder\Raw.mat'], Path1);

        Raw_map = matfile([Files(1:loc) '\processed_data_Folder\Raw.mat'], 'Writable', true);
        clear Y % Yr
        
        % Preallocate matfile handles
        Handles = cell(num_planes, 1);
        
        % Create individual RAw files for each plane
        if num_planes > 1
            
            Y = zeros(X_Axis, Y_Axis, Length_Recording, 'uint16');
            savefast([Files(1:loc) '\processed_data_Folder\Raw_1.mat'], 'Y');  
            Handles{1} = matfile([Files(1:loc) '\processed_data_Folder\Raw_1.mat'], 'Writable', true);
            
            for z = 2:num_planes
                Path1 = [Files(1:loc) '\processed_data_Folder\Raw_' num2str(z) '.mat'];
                copyfile([Files(1:loc) '\processed_data_Folder\Raw_1.mat'], Path1);
                Handles{z} = matfile([Files(1:loc) '\processed_data_Folder\Raw_' num2str(z) '.mat'], 'Writable', true);
            end  
        else 
            copyfile([Files(1:loc) '\processed_data_Folder\Raw.mat'], ...
                [Files(1:loc) '\processed_data_Folder\Raw_1.mat']);
            Handles{1} = matfile([Files(1:loc) '\processed_data_Folder\Raw_1.mat'], 'Writable', true);
        end
        
        clear Y
        
    end
    
    % Create List of indices 
    All_Indices = reshape(1:Length_Recording*num_planes, num_planes, Length_Recording);

    if num_planes > 1
        switch Scanning_mode
            case 'Bidirectional'
                All_Indices(:, 2:2:end) = flipud(All_Indices(:, 2:2:end));
        end
    end
    
    % Now load all tif files in batches and either save them in a memory
    % mapped mat file or a concatenated tif file
    Current_Spot = 1;
    if ~NewStacks
        parpool(num_workers);
        for i = 1:Iterations
            Start_time = tic;

            % Get the Patches 
            Start = 1 + Patch_size*(i-1);

            % Now find the right stop frame number
            if i == Iterations
                Stop = size(All_Indices, 2);
            else
                Stop = Patch_size*(i);
            end

            % create a container file for the Images
            Images = cell(numel(Start:Stop), num_planes);

            % Create List of indices 
            Loop_Indices = All_Indices(:, Start:Stop);

            % Iterate over all Planes
            for iii = 1:num_planes
                Image_temp = cell(numel(Start:Stop), 1);
                Indices_load = Loop_Indices(iii,:);            
                Names = {Filenames(Indices_load).name};

                % Load the files in parallel
                parfor ii = 1:numel(Indices_load)
                    Img = loadtiff([Files_Img Names{ii}]);  
                    Image_temp{ii} = Img;
                end

                % Temporarily write Image into file
                Images(:, iii) = Image_temp;
            end

            if ~MemMap
                for ii = 1:num_planes
                    Concatenated_mov = imresize(cat(3, Images{:,ii}), [X_Axis, Y_Axis], 'bicubic');
                    saveastiff(Concatenated_mov, [Files(1:loc) 'Raw_' num2str(ii) '.tif'], options);
                    clear Concatenated_mov Images Image_temp
                end
            else     
                for ii = 1:num_planes

                    Conc_Imgs = imresize(cat(3, Images{:,ii}), [X_Axis, Y_Axis], 'bicubic');  

                    if num_planes > 1
                        Raw_map.Y(:, :, Current_Spot:Current_Spot + size(Conc_Imgs, 3) - 1, ii) = Conc_Imgs;
                    else
                        Raw_map.Y(:, :, Current_Spot:Current_Spot + size(Conc_Imgs, 3) - 1) = Conc_Imgs;
                    end

                    Handles{ii}.Y(:, :, Current_Spot:Current_Spot + size(Conc_Imgs, 3) - 1) = Conc_Imgs;
                end
                Current_Spot = Current_Spot + size(Conc_Imgs, 3);
            end

            Timing_end = toc(Start_time);        
            fprintf('Finished %.f out of %.f Patches, processing took %.f s \n', i, Iterations, Timing_end)
        end
        delete(gcp)
    else
        if num_planes > 1
            error(['Havent recorded multiplane 2P data, therefore this code' ...
                ' doesnt support loading it yet. Please update this section '...
                ' correspondingly.']);
        end
        
        [~, idx] = sort_nat({Filenames.date}); % Order the indices according to recording time
        for k = 1:size(Filenames, 1)
            Start_time = tic;
            DataLoad = loadtiff([Filenames(idx(k)).folder '\' Filenames(idx(k)).name]);
             for ii = 1
                Conc_Imgs = imresize(DataLoad, [X_Axis, Y_Axis], 'bicubic'); % make square in case it isnt
                
                if num_planes > 1
                    Raw_map.Y(:, :, Current_Spot:Current_Spot + size(Conc_Imgs, 3) - 1, ii) = Conc_Imgs;
                else
                    Raw_map.Y(:, :, Current_Spot:Current_Spot + size(Conc_Imgs, 3) - 1) = Conc_Imgs;
                end
                
                Handles{ii}.Y(:, :, Current_Spot:Current_Spot + size(Conc_Imgs, 3) - 1) = Conc_Imgs;
             end
             Current_Spot = Current_Spot + size(Conc_Imgs, 3);
             Timing_end = toc(Start_time);        
             fprintf('Finished %.f out of %.f Image Stacks, processing took %.f s \n', k, size(Filenames, 1), Timing_end)
        end
    end
end