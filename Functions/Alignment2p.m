%% Alignment of multi day 2P Imaging recordings
% Code written by Julian Hinz @ Luthi lab 2019, 2020

function [Output] = Alignment2p(Files)
        
    Output = [];
    for i = 1:size(Files, 2)
    
        % Current Folder    
        Filepath = Files{i};

        % First create Folder for Alignement
        if ~exist([Filepath '\Alignment\'])
           mkdir([Filepath '\Alignment\'])
        end

        % Automatically detect folders that need to be processed
        files = dir(Filepath);% Get a list of all files and folders in this folder.
        dirFlags = [files.isdir]; % Get a logical vector that tells which is a directory.
        subFolders = files(dirFlags); % Extract only those that are directories.

        % Sort the folders and select the correct ones based on length and name

        index_Ses = find(cellfun(@(x) any(strfind(x, 'Session')), {subFolders(:).name}) & ...
            cell2mat(cellfun(@numel, {subFolders(:).name},'UniformOutput',false)) < 20);
        subFolders = subFolders(index_Ses); 
        subFolders = sort_nat({subFolders(:).name});

        % Preallocate variables
        Alignment = struct;
        mult_planes = false;
        Max_session = nan(size(subFolders, 1), 2);


        % Run through all Sessions, load the 
        for ii = 1:size(subFolders, 2)
            try
                load([Filepath '\' subFolders{ii} '\processed_data_Folder\PospP_CNMF.mat'])
                load([Filepath '\' subFolders{ii} '\processed_data_Folder\Corr_Imgs.mat'])
            catch
                disp('Session could not be loaded');
                continue
            end

            temporary_names = fieldnames(neurons_f);
            Corr_Imgs = cell(size(temporary_names, 1), 1);
            Components = cell(size(temporary_names, 1), 1);

            for z = 1:size(temporary_names, 1)
                Corr_Imgs{z} = obj.(['Cor_Images_' num2str(z)]).Corr_Img(:,:,1); 
                Components{z, 1} = neurons_f.(temporary_names{z}).A;
                Components{z, 2} = neurons_f.(temporary_names{z}).Centroids;
            end

            Neuron_mult.(['Session_' num2str(ii)]).neuron = neurons_f;
            Neuron_mult.(['Session_' num2str(ii)]).Spatial = obj;

            Alignment.(['Session_' num2str(ii)]).Correlation_Imgs = cat(3, Corr_Imgs{:});
            Alignment.(['Session_' num2str(ii)]).Components = Components;

            Max_session(ii) = size(temporary_names, 1);

            if size(temporary_names, 1) ~= 1
               mult_planes = true;
            end

        end

        % If there is at least one mult sessions, we need to make sure that we
        % try to align the rights planes across sessions

        Alignment.mult_planes = mult_planes;
        Alignment.Max_planes = Max_session(~isnan(Max_session));

        % Now save this, so we can have fast access later
        savefast([Filepath '\Alignment\SpatialInfo_collected.mat'], 'Alignment');
        savefast([Filepath '\Alignment\Sessions_combined.mat'], 'Neuron_mult');    
    end


    % First we need to make sure that (1) we only try to align Sessions that we
    % actually want to align and (2) that we align the right planes together

    temp = [];

    for i = 1:size(Files, 2)

        % Current Folder    
        Filepath = Files{i};

        % Load Alignment Info
        load([Filepath '\Alignment\SpatialInfo_collected.mat']);

        % Offer to concatenate all sessions, if all Sessions have the same
        % numbers of Planes imaged otherwise go either into displaying the
        % indv. Sessions

        if range(Alignment.Max_planes) == 0
            fprintf('Do you want to concatenate all sessions ? - Press Enter for yes, press any key and then enter for no:    ');
            temp = input('', 's');
        end 

        % If only one plane was imaged and all sessions should be concatenated,
        % terminate loop here

        if isempty(temp) && ~mult_planes
           Align_Ses = logical(ones(size(Alignment.Max_planes, 2), 1)); 
           Alignment.Align_Ses = Align_Ses;
           Alignment.Order = mat2cell(ones(size(Alignment.Max_planes, 2), 1), ones(size(Alignment.Max_planes, 2), 1), 1);
           continue
        end

        % Now check for which sessions should be kicked out and also, reorder
        % if the planes are not in the right order

        tmp_names = fieldnames(Alignment);
        index_Ses = find(cellfun(@(x) any(strfind(x, 'Session')), tmp_names));
        tmp_names = tmp_names(index_Ses); 

        planes_max = max(Alignment.Max_planes);
        Align_Ses = zeros(size(tmp_names, 1), 1);
        Order_coll = cell(size(tmp_names, 1), 1);

        % Create the two figures
        f = figure;
        g = figure;

        % Run through all the sessions and check whether they should be aligned
        for z = 1:size(tmp_names, 1)   

            % Load Session for easier access
            tmp_Aligninfo = Alignment.(tmp_names{z});
            Ses_del = false;

            % Now plot and replot the figure until the optimal Order is reached
            Run_Time = 0;
            tmp_Order = [1 2 3 4];
            approval = false;

            while ~Ses_del && ~approval

                % Make the loop flexible for multiple runs
                Run_Time = Run_Time + 1;

                % Plot the correlation images to check before alignment, also
                % link axis to enable synchronous zoom
                figure(g)
                if ~any(Align_Ses(1:z))
                    for zz = 1:size(tmp_Aligninfo.Correlation_Imgs, 3)
                        Handles{zz} = subplot(size(tmp_Aligninfo.Correlation_Imgs, 3), 1, zz);
                        imagesc(tmp_Aligninfo.Correlation_Imgs(:, :, tmp_Order(zz))); axis off; 
                        colormap bone
                    end

                    linkaxes([Handles{:}])

                else

                    idx_nonzero = find(Align_Ses(1:z) ~= 0);
                    indices = [1:1:planes_max];

                    for zz = 1:size(tmp_Aligninfo_last.Correlation_Imgs, 3)
                        Handles2{zz} = subplot(planes_max, 2, indices(zz)*2 - 1);
                        imagesc(tmp_Aligninfo_last.Correlation_Imgs(:, :, Order_coll{idx_nonzero(end)}(zz))); axis off; 
                        colormap bone
                    end

                    for zz = 1:size(tmp_Aligninfo.Correlation_Imgs, 3)
                        Handles{zz} = subplot(planes_max, 2, indices(zz)*2);
                        imagesc(tmp_Aligninfo.Correlation_Imgs(:, :, tmp_Order(zz))); axis off; 
                        colormap bone
                    end

                    linkaxes([Handles{:} Handles2{:}]);

                end

                % Now, Check which order and whether Sessions should be
                % processed
                if Run_Time == 1
                    fprintf(['Do you want to delete this session from concatenation procedure ?' ...
                        ' - Press Enter for No, press any key and then enter for yes:      ']);
                    temp = input('', 's');

                    if ~isempty(temp)
                        Ses_del = true;
                        continue
                    end

                    right_numbers = false;

                    while ~right_numbers
                        fprintf('Please specify the right Order of Planes with numbers seperated by spaces: ');
                        tmp_input = input('', 's');
                        tmp_Order = str2num(tmp_input);

                        if ~isempty(tmp_Order) && isnumeric(tmp_Order) && ...
                                size(tmp_Order, 2) == size(tmp_Aligninfo.Correlation_Imgs, 3)
                            right_numbers = true;
                        end
                    end
                else
                    fprintf('Are you happy with the order ? - Press Enter for Yes, press any key and then enter for No:');
                    temp = input('', 's');

                    if isempty(temp) && ~Ses_del
                        Align_Ses(z) = 1;
                        approval = true;
                    elseif Ses_del
                        approval = true;                    
                    else
                        Run_Time = 0;
                    end
                end

            end

            disp('############################################################')
            disp('############################################################')

            % Plot the results in the overview figure
            if ~Ses_del
                figure(f)
                numbers = [z z+size(tmp_names, 1) ...
                        z+2*size(tmp_names, 1) z+3*size(tmp_names, 1)];
                for zz = 1:size(tmp_Aligninfo.Correlation_Imgs, 3)
                    subplot(planes_max, size(tmp_names, 1), numbers(zz)); % (z-1)*planes_max + zz
                    imagesc(tmp_Aligninfo.Correlation_Imgs(:, :, tmp_Order(zz))); axis off; 
                    colormap bone
                end

                Order_coll{z} = tmp_Order;
                tmp_Aligninfo_last = tmp_Aligninfo;
            end

            % Remove the Axes for the next iteration
            try
                cellfun(@cla, Handles);
                clear Handles

                cellfun(@cla, Handles2);
                clear Handles2 
            end

        end

        % Save the results
        saveas(f, [Filepath '\Alignment\Order of planes.png'])
        saveas(f, [Filepath '\Alignment\Order of planes.fig'])

        close(g)
        close(f)

        Alignment.Align_Ses = logical(Align_Ses);
        Alignment.Order = Order_coll;

        savefast([Filepath '\Alignment\SpatialInfo_collected.mat'], 'Alignment')
        savefast([Filepath '\Alignment\Order_Planes.mat'], 'Order_coll');
        savefast([Filepath '\Alignment\Sessions_tobealigned.mat'], 'Align_Ses');
    end


    % After collecting the Information about which sessions should be
    % processed, now enter the actual Alignment stage

    %%%%%%%  Registration over days
    reference_Session = 1; % Reference Session for registration
    max_Distance = 7; % in pix, distance that is maximally allowed for still belonging to the same group
    analysis_Type = 'pairwise'; % 'clustering'

    [optimizer, metric] = imregconfig('multimodal');
    optimizer.MaximumIterations = 1000;
    optimizer.InitialRadius = 0.00001;

    for i = 1:size(Files, 2)

        % Current Folder    
        Filepath = Files{i};

        % Load Alignment Info
        load([Filepath '\Alignment\SpatialInfo_collected.mat']);

        Ses_toalign = find(Alignment.Align_Ses == 1);
        Order_Align = cell2mat(Alignment.Order(Ses_toalign));
        Image_Size = size(Alignment.(['Session_' num2str(Ses_toalign(reference_Session))]).Correlation_Imgs);
        
        %%% Ask whether user wants to specify Region for Registration
        fprintf(['Do you want to do manual or automatic Registration ? -' ... 
            'Press Enter for manual, press any key and then enter for No:']);
        temp = input('', 's');
        
        if isempty(temp)
            Registration_Manual = true;
        else
            Registration_Manual = false;
        end
        
        if ~Registration_Manual
            %%% Ask whether user wants to specify Region for Registration
            fprintf(['Do you want to specify Region for Registration ? -' ... 
                'Press Enter for Yes, press any key and then enter for No:']);
            temp = input('', 's');

            Selection = cell(size(Order_Align, 2), 1); 
            if isempty(temp)
               for zz = 1:size(Order_Align, 2)
                   [Selection{zz}] = Define_ROIS_2P(Alignment.(['Session_' num2str(Ses_toalign(reference_Session))]). ...
                    Correlation_Imgs(:, :, Order_Align(1, zz)), 1);
               end
            else
                for zz = 1:size(Order_Align, 2)
                    tmp_IMGS = Alignment.(['Session_' num2str(Ses_toalign(reference_Session))]). ...
                        Correlation_Imgs(:, :, Order_Align(1, zz));
                    Selection{zz} = [{1:size(tmp_IMGS, 1)}, {1:size(tmp_IMGS, 2)}];  
                end
            end
        end
        
        % Run through the planes
        for zz = 1:size(Order_Align, 2)
            
            Corr_Imgs = nan(Image_Size(1), Image_Size(2), numel(Ses_toalign));
            Corr_Imgs(:, :, 1) = Alignment.(['Session_' num2str(Ses_toalign(reference_Session))]). ...
                Correlation_Imgs(:, :, Order_Align(1, zz));
            reg = cell(numel(Ses_toalign) - 1, 1);
            
            
            if ~Registration_Manual
                Extract_Selection = Selection{zz};

                bw_x = Extract_Selection{1, 1};
                bw_y = Extract_Selection{1, 2};

                % First use the Correlation Images to calculate the rigid body
                % shift between sessions related to the session immediatly prior,
                % but all in reference to the first

                for zzz = 2:numel(Ses_toalign) % add the ability to align sessions to a random session
                    tmp_Aligninfo = Alignment.(['Session_' num2str(Ses_toalign(zzz))]);
                    Temp_Corr_Imgs = tmp_Aligninfo.Correlation_Imgs(:, :, Order_Align(zzz, zz));
                    
                    reg{zzz-1} = imregtform(Temp_Corr_Imgs(bw_y, bw_x), Corr_Imgs(bw_y, bw_x, ...
                        zzz-1), 'rigid', optimizer, metric);
                    Corr_Imgs(:,:,zzz) = imwarp(Temp_Corr_Imgs, reg{zzz-1}, ...
                        'OutputView', imref2d( size(Temp_Corr_Imgs) ));
                end
            else
                
                Reference = imref2d(size(Corr_Imgs(:,:,1)));    
                
                for zzz = 2:numel(Ses_toalign)
                    tmp_Aligninfo = Alignment.(['Session_' num2str(Ses_toalign(zzz))]);
                    Temp_Corr_Imgs = tmp_Aligninfo.Correlation_Imgs(:, :, Order_Align(zzz, zz));
                    
                    [Ref_Points_moving, Ref_Points_fixed] = cpselect(Temp_Corr_Imgs,...
                        Corr_Imgs(:,:,zzz-1), 'Wait',true);
                    
                    reg{zzz - 1} = fitgeotrans(Ref_Points_moving, Ref_Points_fixed, ...
                        'nonreflectivesimilarity');
                    
                    %relate intrinsic and world coordinates
                    Corr_Imgs(:, :, zzz) = imwarp(Temp_Corr_Imgs, reg{zzz - 1}, ...
                        'OutputView', Reference);
                    
                end
                
            end
            
            %%%%% Same Process in both cases
            
            % Now transform the neural footprints according to the calculated
            % transformation

            Coll_Spatial = cell(numel(Ses_toalign), 1);

            for zzz = 1:numel(Ses_toalign)
                Comp_temp = Alignment.(['Session_' num2str(Ses_toalign(zzz))]). ...
                    Components{Order_Align(zzz, zz), 1};
                Temporary_befTransformation = reshape(Comp_temp, Image_Size(1), Image_Size(2), ...
                    size(Comp_temp, 2));
                After_Transformation = nan(size(Temporary_befTransformation));
                
                Correct = true(size(Temporary_befTransformation, 3), 1);
                
                if zzz == 1 
                    After_Transformation = Temporary_befTransformation;
                else
                    for zzzz = 1:size(Temporary_befTransformation, 3)
                        tmp = imwarp(Temporary_befTransformation(:, :, zzzz), ...
                            reg{zzz-1}, 'OutputView', imref2d( size(Temporary_befTransformation(:,:,1)) ));
                        if ~any(any(tmp))
                            Correct(zzzz) = false;
                            continue
                        else
                            After_Transformation(:,:,zzzz) = tmp;
                        end
                        
                    end 
                end
                
                Coll_Spatial{zzz} = After_Transformation(:, :, Correct);
                savefast([Filepath '\Alignment\ComponentsinAlignment_Plane' num2str(zz) '.mat'], 'Correct')
            end

            % Save all the Alignment results
            savefast([Filepath '\Alignment\Registration_across_Ses_transformation_Plane' num2str(zz) '.mat'], 'reg')
            savefast([Filepath '\Alignment\Aligned_CorrImgs_Plane' num2str(zz) '.mat'], 'Corr_Imgs')
            savefast([Filepath '\Alignment\Transformed_SpatialComponents_Plane' num2str(zz) '.mat'], 'Coll_Spatial')

            % Now, match the spatial components across the the Sessions
            [OutStruct] = matchObjBtwnTrials2(Coll_Spatial, max_Distance, ...
                reference_Session, analysis_Type);

            %%% Post Processing

            try
                Num_Ses_aligned = cumsum(transpose(OutStruct.globalIDs~=0))';
            catch
                disp('No automatic Alignment possible !')
                continue
            end

            b = hist(Num_Ses_aligned(:,end), [1:size(Num_Ses_aligned, 2)]);
            g = figure;
            bar(b, 'k'); box off
            xlabel('Number of aligned Sessions'); ylabel('# Cells')
            saveas(gcf, [Filepath '\Alignment\Alignment_Plane_' num2str(zz) '.png'])
            close(g)

            [Ses, index] = sortrows(Num_Ses_aligned, 'descend');
            All_Ses = max(find(Ses(:,end) == size(Ses, 2)));

            savefast([Filepath '\Alignment\Component_Alignment_Plane' num2str(zz) '.mat'], 'OutStruct');
            close all
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Old Code that was used previously and is now integrated in larger
%%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test = nan(512, 512, numel(3:21));
% 
% for zz = 3:21
%     tmp_Aligninfo = Alignment.(tmp_names{zz});
%     Test(:,:, zz-2) = tmp_Aligninfo.Correlation_Imgs(:, :, 1);
% end
% 
% % [Position_all] = Define_ROIS_mini(Test(:, :, 1), 1);
% 
% Test2 = nan(512, 512, size(Test, 3));
% 
% [optimizer, metric] = imregconfig('multimodal');
% optimizer.MaximumIterations = 200;
% 
% Test2(:,:,1) = Test(:,:,1);
% reg = cell(size(Test, 3) - 1, 1);
% 
% for z = 2:size(Test, 3)
%     [Test2(:,:,z), ~] = imregister(Test(:,:,z), Test2(:,:, z-1), 'rigid', optimizer, metric);
%     reg{z-1} = imregtform(Test(:,:,z), Test2(:,:, z-1), 'rigid', optimizer, metric);
% end
% 
% Coll_Spatial = cell(numel(Ses_toalign)-1, 1);
% 
% for zzz = 2:numel(Ses_toalign)
%     
%     Comp_temp = Alignment.(['Session_' num2str(Ses_toalign(zzz))]).Components{Order_Align(zz), 1};
%     Temporary_befTransformation = reshape(Comp_temp, Image_Size(1), Image_Size(2), ...
%         size(Comp_temp, 2));
%     After_Transformation = nan(size(Temporary_befTransformation));
%     
%     if zzz == 2 % 1 
%         After_Transformation = Temporary_befTransformation;
%     else
%         for zzzz = 1:size(Temporary_befTransformation, 3)
%             After_Transformation(:,:,zzzz) = imwarp(Temporary_befTransformation(:,:,zzzz), ...
%                 reg{zzz-2}, 'OutputView', imref2d( size(Temporary_befTransformation(:,:,1)) ));
%         end 
%     end
%     Coll_Spatial{zzz-1} = After_Transformation;
% end
% 
% 
% 
% 
% % for z = 2:18
% %     Test2(:,:,z) = imregister(Test(:,:,z), Test(:,:, 16), 'rigid', optimizer, metric);
% % end
% 
% 
% 
% writeTiff(Corr_Imgs, '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Matregisternew.tif');











