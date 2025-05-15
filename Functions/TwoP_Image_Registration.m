%% Analyze 2P Imaging

function [F, Templates_rg, Shifts_All] = ...
    TwoP_Image_Registration(Image_location, MemMap, T_DS_Factor, S_DS_Factor,...
    DS_Factor, num_planes, Selection, Patch_size, Run_time, num_workers)
        
    Shifts_All = cell(numel(num_planes), 1);
    Templates_rg = cell(numel(num_planes), 1);
    
    Single_Plane = numel(num_planes) == 1;
    
    % Start the parallel pool
    parpool(num_workers);
    
    for zz = 1:numel(num_planes)

        % Set parameters
        if MemMap 
            ext = '.mat';
            % Set parameters for Motion correction
            Raw_mat = matfile([Image_location '\processed_data_Folder\Raw_' num2str(num_planes(zz)) ext]);
            Size_sm = size(Raw_mat.Y);
            options_rg = NoRMCorreSetParms('d1', Size_sm(1), 'd2', Size_sm(2), ...
                'bin_width',150,'max_shift',15,'us_fac',50, 'output_type', 'memmap', ...
                'use_windowing', true, 'phase_flag', true, 'mem_batch_size', 5000);  
        else          
            error('Script does not support tiff anymore !')
        end

        % NormCorre registration
        [~, shifts_g, Templates_tmp, ~, ~] = normcorre_batch_adjusted( ...
            [Image_location '\processed_data_Folder\Raw_' num2str(num_planes(zz)) ext],...
            options_rg, [], true, zz, Single_Plane, Selection);
        
        % Shifts applied
        shifts_r = squeeze(cat(3, shifts_g(:).shifts));
        
        % If the only part of the frame was used for MC, apply the shifts
        % to the full movie
        if ~isempty(Selection)
            
            Raw_mat = matfile([Image_location '\processed_data_Folder\Raw_' num2str(num_planes(zz)) ext]);
            MC_mat = matfile([Image_location '\processed_data_Folder\MC_r' ext], 'Writable', true);
            
            Iterations = ceil(size(shifts_r, 1)/Patch_size);
            
            for k = 1:Iterations
                
                % Get the Patches 
                Start = 1 + Patch_size*(k-1);

                % Now find the right stop frame number
                if k == Iterations
                    Stop = size(shifts_r, 1);
                else
                    Stop = Patch_size*(k);
                end
                
                % Load Patch and create empty container, then loop over the
                % data and apply shifts
                Y_tmp = Raw_mat.Y(:,:,Start:Stop);
                
                class = whos('Y_tmp');
                Y_tmp_reg = zeros(size(Y_tmp), class.class);
                shifts_tmp = shifts_r(Start:Stop, :);
                
                parfor kk = 1:size(Y_tmp, 3)
                    Y_tmp_reg(:,:,kk) = imtranslate(Y_tmp(:,:,kk), [shifts_tmp(kk, 2), shifts_tmp(kk, 1)]);
                end
               
                if Single_Plane
                    MC_mat.Y(:, :, Start:Stop) = Y_tmp_reg;
                else
                    MC_mat.Y(:, :, Start:Stop, zz) = Y_tmp_reg;
                end
                                
                clear Y_tmp Y_tmp_reg
                
            end
        end
        
        Shifts_All{zz} = [shifts_r sum(abs(shifts_r), 2)];
        clear shifts_g
        
        Templates_rg{zz} = Templates_tmp;
        
        fprintf('\nFinished Plane %d out of %d \n', zz, numel(num_planes));
    end
    
    delete(gcp);
        
    % First, Save the DS video 
    Y = [];
    K = matfile([Image_location '\processed_data_Folder\MC_r' ext]);
    Sizes = size(K.Y);

    DS_overhang = mod(Sizes(3), T_DS_Factor);

    if Single_Plane

        if T_DS_Factor > 1 || S_DS_Factor < 1

            % Do Temporal and Spatial DS if requested
            DS_Dat = reshape(squeeze(mean(reshape(K.Y(:, :, 1:end-DS_overhang),...
                [Sizes(1)*Sizes(2), T_DS_Factor, (Sizes(3)-DS_overhang)/T_DS_Factor]), 2, 'native')), ...
                [Sizes(1), Sizes(2), (Sizes(3)-DS_overhang)/T_DS_Factor]);
            Sizes_n = size(DS_Dat);
            DS_Dat = imresize(DS_Dat, [Sizes_n(1)*S_DS_Factor, Sizes_n(2)*S_DS_Factor], 'bicubic');
            
            if Run_time == 1
                % Now change the Timing File accroding to Temporal DS
                Timing_Info = struct;
                load([Image_location '\processed_data_Folder\Timing_Information.mat'])
                Timing_Info.FS = Timing_Info.FS/T_DS_Factor;
                if Sizes(3) == size(Timing_Info.Timing, 1)
                    Size_Time = size(Timing_Info.Timing(1:end-DS_overhang, :));
                    Timing_Info.Timing = squeeze(mean(reshape(Timing_Info.Timing(1:end-DS_overhang, :), T_DS_Factor, Size_Time(1)/T_DS_Factor, ...
                        Size_Time(2)), 1));
                else
                    Size_Time = size(Timing_Info.Timing(1:end, :));
                    Timing_Info.Timing = squeeze(mean(reshape(Timing_Info.Timing(1:end, :), T_DS_Factor, Size_Time(1)/T_DS_Factor, ...
                        Size_Time(2)), 1));
                end
                savefast([Image_location '\processed_data_Folder\Timing_Information.mat'], 'Timing_Info')
            end
        else
            DS_Dat = K.Y;        
        end

        savefast([Image_location '\processed_data_Folder\DS_Dat' ext], 'DS_Dat')
        clear DS_Dat

        DS_overhang = mod(Sizes(3), DS_Factor);
        DS_Dat_vis = uint16(reshape(squeeze(mean(reshape(K.Y(:, :, 1:end-DS_overhang), [Sizes(1)*Sizes(2), DS_Factor, (Sizes(3) - DS_overhang)/DS_Factor]), 2)), ...
            [Sizes(1), Sizes(2), (Sizes(3) - DS_overhang)/DS_Factor]));

        Max_val = max(DS_Dat_vis(:));
        savefast([Image_location '\processed_data_Folder\DS_Dat_vis' ext], 'DS_Dat_vis')
        clear Y_temp

        F(numel(1:size(DS_Dat_vis, 3))) = struct('cdata',[],'colormap',[]);
        fig = figure('position', [100, 50, 1300, 1300]); 
        for z = 1:size(DS_Dat_vis, 3)
            imshow(DS_Dat_vis(:,:,z),[0 Max_val*0.7]);                
            drawnow; 
            F(z) = getframe(fig);
        end 

        close(fig)

    else  

        if T_DS_Factor > 1 || S_DS_Factor < 1

            % Do Temporal and Spatial DS if requested
            DS_Dat = reshape(squeeze(mean(reshape(K.Y(:, :, 1:end-DS_overhang, :), [Sizes(1)*Sizes(2), T_DS_Factor, (Sizes(3)-DS_overhang)/T_DS_Factor, Sizes(4)]), 2, 'native')), ...
                [Sizes(1), Sizes(2), (Sizes(3)-DS_overhang)/T_DS_Factor, Sizes(4)]);

            Sizes_n = size(DS_Dat);
            DS_Dat = imresize(DS_Dat, [Sizes_n(1)*S_DS_Factor, Sizes_n(2)*S_DS_Factor], 'bicubic');
            
            if Run_time == 1
                % Now change the Timing File accroding to Temporal DS
                Timing_Info = struct;
                load([Image_location '\processed_data_Folder\Timing_Information.mat'])
                Timing_Info.FS = Timing_Info.FS/T_DS_Factor;
                Size_Time = size(Timing_Info.Timing(1:end-DS_overhang, :, :));
                Timing_Info.Timing = squeeze(mean(reshape(Timing_Info.Timing(1:end-DS_overhang, :, :), T_DS_Factor, Size_Time(1)/T_DS_Factor, ...
                    Size_Time(2), Size_Time(3)), 1));

                savefast([Image_location '\processed_data_Folder\Timing_Information.mat'])
            end
            
        else
            DS_Dat = K.Y;        
        end

        savefast([Image_location '\processed_data_Folder\DS_Dat' ext], 'DS_Dat')
        clear DS_Dat

        DS_overhang = mod(Sizes(3), DS_Factor);
        DS_Dat_vis = uint16(reshape(squeeze(mean(reshape(K.Y(:, :, 1:end-DS_overhang, :), [Sizes(1)*Sizes(2), DS_Factor, (Sizes(3)-DS_overhang)/DS_Factor, Sizes(4)]), 2)), ...
            [Sizes(1), Sizes(2), (Sizes(3)-DS_overhang)/DS_Factor, Sizes(4)]));

        Max_val = max(DS_Dat_vis(:));
        savefast([Image_location '\processed_data_Folder\DS_Dat_vis' ext], 'DS_Dat_vis')
        clear Y_temp

        F(numel(1:size(DS_Dat_vis, 3))) = struct('cdata',[],'colormap',[]);
        fig = figure('position', [100, 50, 1300, 1300]); 
        for z = 1:size(DS_Dat_vis, 3)
            for zz = 1:numel(num_planes)
                subplot(ceil(numel(num_planes)/3), 3, zz)
                imshow(DS_Dat_vis(:,:, z, zz),[0 Max_val*0.7]);                
            end 
            drawnow; 
            F(z) = getframe(fig);
        end 

        close(fig)
    end

end
 
%% Old code, based on loading all in memory



%     MC_Movies = cell(size(Filenames, 1), 2);

%         Data_Temp =  loadtiff([Image_location '\' Filenames(i).name]);
%         MRDFT = imresize(Data_Temp(:,:,1:T_DS:end), [size(Data_Temp, 1)*S_DS size(Data_Temp, 2)*S_DS], ...
%              'bicubic');
%          
%         [~, ~, T] = size(MRDFT);    % dimensions of file
%         MRDFT = MRDFT - min(MRDFT(:));      % remove negative offset
% 
%         minY = quantile(MRDFT(1:1e7),0.0005);
%         maxY = quantile(MRDFT(1:1e7),1-0.0005);
%         
%         if vis
%             %  view data
%             figure;play_movie({MRDFT},{'raw data'},minY,maxY);
%         end
        
        % [a, b, Indicesll] = Define_ROIS(MRDFT(:,:,1), []);
        

%         Shifts_summarized = [];

%         [M_rg, shifts_rg,template_rg] = normcorre_batch(MRDFT, options_rg);        
%         MC_Movies{i, 1} = M_rg;
%         

%         if vis
%             view data
%             tsub = 5;   % downsampling factor (only for display purposes)
%             Y_sub = downsample_data(MRDFT,'time',tsub);
%             M_rgs = downsample_data(M_rg,'time',tsub);
% 
%             play_movie({Y_sub,M_rgs},{'raw data','rigid'},minY,maxY);
%         end
        
        % perform non-rigid motion correction    
        % parameters motion correction
        
        % 'd1','d2': size FOV movie
        % 'grid_size','overlap_pre': parameters regulating size of patch (size patch ~ (grid_size + 2*overlap_pre))
        % 'mot_uf': upsampling factor of the grid for shift application
        % 'bin_width': how often to update the template
        % 'max_shift': maximum allowed rigid shift
        % 'max_dev': maximum deviation allowed for each patch from the rigid shift value


        




%% Seems to be worthless, Ill try with full CNMF protocol    
% Tried to identify Images with shifts ... I just recorded better data
%         Transformation = nan(1, size(M_rg, 3));
%         Transformation2 = nan(1, size(M_rg, 3));

%         parfor kkkk = 2:size(M_rg, 3)
%             Temp = normxcorr2(template_rg(a{1,2}, a{1,1}), M_rg(a{1,2}, a{1,1},kkkk));
%             Transformation(:,kkkk) = max(abs(Temp(:)));
%             Temp2 = normxcorr2(M_rg(a{1,2}, a{1,1},kkkk - 1), M_rg(a{1,2}, a{1,1},kkkk));
%             Transformation2(:,kkkk) = max(abs(Temp2(:)));
% %             Transformed_imgs(:,:,kkkk) = imwarp(M_rg(:,:,kkkk), Temp, 'OutputView', Rfixed); 
%         M_rg_gray = mat2gray(M_rg);implay(M_rg_gray);
%         figure;plot(Transformation)
%         figure;plot((Transformation2))
%         
%         end    
    
%     % Check why motion correction not working.. Probably something with the
%     % Imgaussfilter, also need to check how to identify Images that are 
%     % affected by motion, might be possible by looking at amount of shift 
%     for i = 1:size(Filenames, 1)
%         
%         Data_Temp =  loadtiff([Image_location '\' Filenames(i).name]);
%         Data{i} = imresize(Data_Temp(:,:,1:T_DS:end), [size(Data_Temp, 1)*S_DS size(Data_Temp, 2)*S_DS], ...
%              'bicubic');
%          
%         disp('Finished loading movie !')
%         
%         Position_all = MCR_Areas{i};
%         
%         for ii = 1:size(Position_all, 1)
% 
%             % Do Motion Correction, therefore first create Template
%             if ii == 1
%                 MRDFT_Template = Data{i}(:,:,1:100);
%             else
%                 MRDFT_Template = MRDFT(:,:,1:100);
%             end
% 
%             Start_time = tic; count = 0; Shifts = 1;
%             while mean(abs(Shifts(:))) > Tolerance_motion && count < max_reps_motioncorr                 
%                 count = count + 1;
%                 [MRDFT_Template, Shifts] = Image_Registration_2P_rigid(MRDFT_Template, MRDFT_Template(:,:,1), 2, S_DS, Position_all{ii,1}, Position_all{ii,2});                
%             end
%             Time_final = toc(Start_time);
%             fprintf('Template correction converged after %d repetitions and took %e s\n', count, Time_final);
% 
%             Template = median(MRDFT_Template, 3); % create new template
%             clear MRDFT_Template
% 
%             % Then run through the Data and do motion correction
% 
%             if ii == 1
%                 MRDFT = Data{i};
%             end
% 
%             Start_time = tic; count = 0; Shifts = 1;
%             Temp_Buf_Results = cell(1);
% 
%             while mean(abs(Shifts(:))) > Tolerance_motion && count < max_reps_motioncorr                 
% 
%                 count = count + 1;
%                 [MRDFT, Shifts] = Image_Registration_2P_rigid(MRDFT, Template, 2, S_DS, Position_all{i,1}, Position_all{i,2});
% 
%                 if mean(mean(abs(Shifts - mean(Shifts, 1)))) < Tolerance_motion && count > 1
%                     MRDFT = Temp_Buf_Results{1};
%                     break
%                 else
%                     Temp_Buf_Results{1} = MRDFT;
%                 end
% 
%             end
%             Time_final = toc(Start_time);
%             fprintf('Motion correction converged after %d repetitions and took %e s\n', count, Time_final);
%         end
%         
%         MC_Movies{i} = MRDFT;
%         Templates{i} = Template;
%         
%     end
%     
%     
%     figure;
%     for i = 1:5:size(MRDFT, 3)
%         imagesc((MRDFT(:,:,i)));box off;set(gca, 'Xtick', [], 'Ytick', []); drawnow
%     end
%     Cn = correlation_image(mat2gray(MRDFT),8,512,512, 0);
%     
%     
%     figure;
%     for i = 1:size(MRDFT_Template, 3)
%         imagesc(MRDFT_Template(:,:,i));box off;set(gca, 'Xtick', [], 'Ytick', []); drawnow
%     end
%    [~, ~, Indicesll] = Define_ROIS(MRDFT(:,:,1), []);
%         
%     lin_MRDFT = reshape(imgaussfilt(MRDFT, 3), [size(MRDFT, 1)*size(MRDFT, 2) size(MRDFT, 3)]);
%     All_corrs = corrcoef(lin_MRDFT(Indicesll{1, 1}, :));
%     nmse = nan(size(MRDFT, 3), 1);
%     Template_lin = reshape(Template, 512*512);
%     
%     for k = 1:size(MRDFT, 3)
%        nmse(k) = norm(double(imabsdiff(Template_lin(Indicesll{1, 1}), MRDFT(Indicesll{1, 1}))),'fro');
%     end
%     
%     Radius = 3;
%     Correlations = nan(Radius*2+1, size(MRDFT, 3) - 2*Radius);
%     Correlations2 = nan(2, size(MRDFT, 3) - 2*Radius);
%     for k = 1+Radius:size(MRDFT, 3)-Radius
%         Correlations2(:, k-Radius) = [mean(All_corrs(k, k-Radius)) mean(All_corrs(k, k+Radius))];        
%     end
%     figure;imagesc(Correlations2', [0.7 1]);box off;set(gca, 'Xtick', [], 'Ytick', []);
%     
%     Indices_Low_corr = find(Correlations(5, :) < 0.9 & Correlations(7, :) < 0.9) + Radius;
%     
%     for i = 1:size(Indices_Low_corr, 2)
%         figure;
%         subplot(1,3,1);imagesc(MRDFT(:,:,Indices_Low_corr(i)-1));box off;set(gca, 'Xtick', [], 'Ytick', []); drawnow
%         subplot(1,3,2);imagesc(MRDFT(:,:,Indices_Low_corr(i)));box off;set(gca, 'Xtick', [], 'Ytick', []); drawnow
%         subplot(1,3,3);imagesc(MRDFT(:,:,Indices_Low_corr(i)+1));box off;set(gca, 'Xtick', [], 'Ytick', []); drawnow
%     end
%       
%     cleaned = MRDFT(:,:,Correlations(5, :) < 0.93 & Correlations(7, :) < 0.93);
%     
%     figure;
%     for i = 1:2:size(cleaned, 3)
%         imagesc((cleaned(:,:,i)));box off;set(gca, 'Xtick', [], 'Ytick', []); drawnow
%     end
%     
%     
%     Average_Fluorescence = mean(reshape(MRDFT, [size(MRDFT, 1)*size(MRDFT, 2) size(MRDFT, 3)]), 1);
%     
%     
% %     
% %     [Registered_Images, Shifts] = Image_Registration_2P_rigid(Data{3}, Data{3}(:,:,1), 2, 190:250);
% % 
% % 
% %     options_rigid = NoRMCorreSetParms('d1',size(Data{1},1),'d2',size(Data{1},2),'bin_width',50,'max_shift',40,'us_fac',50, 'shifts_method', 'FFT');
% %     [M_final,shifts,template,options,col_shift] = normcorre(Data{i},options_rigid, Data{i}(:,:,1));



