%% Manualy update CNMF components and then run the final part of the Pipeline

function CNMF_refinement(Filepath, RegionSub, manual_Selection)

    % Load CNMF results   
    if exist([Filepath '\PospP_CNMF.mat']) ~= 0
        disp('Folder was already processed, aborted !')
        return
    elseif exist([Filepath '\Results_CNMF.mat']) ~= 0
        load([Filepath '\Results_CNMF.mat']);
    else
        disp('Folder was not processed, aborted !')
        return
    end
    
    % Load Timing file
    load([Filepath '\Timing_Information.mat'])
    FS = Timing_Info.FS;
    
    % Also load the the correlation file and pass it to the manual
    % selection
    load([Filepath '\Corr_Imgs.mat']);
    
    fields = fieldnames(neuron);
    
    Refinement = cell(size(fields, 1), 2);
    for z = 1:size(fields, 1)
        [delete, del_corr, Data] = manual_Selection_CNMF(neuron.(fields{z}), FS, ...
            obj.(['Cor_Images_' num2str(z)]).Corr_Img(:,:,4), RegionSub); 
        Refinement{z, 1} = delete;
        Refinement{z, 2} = del_corr;
        neurons.(['Plane_' num2str(z)]) = Data;
    end
    
    savefast([Filepath '\Flagged_Components.mat'], 'Refinement');
    savefast([Filepath '\Intermediate_CNMF.mat'], 'neurons');
    
    for z = 1:size(fields, 1)
        
        % Load Current Plane
        tmp_dat = neurons.(['Plane_' num2str(z)]);
        
        % First, delete sorted components
        [tmp_neuron] = update_CNMF_Results(tmp_dat, Refinement{z, 2});

        % detrend fluorescence and extract DF/F values

        tmp_neuron.options.df_window = 1000; 
        [F_dff,F0] = detrend_df_f(tmp_neuron.A,tmp_neuron.b,tmp_neuron.C,tmp_neuron.f, ...
            tmp_neuron.YrA,tmp_neuron.options);

        % deconvolve data

        nNeurons = size(F_dff,1);
        C_dec = zeros(size(F_dff));
        S = zeros(size(F_dff));
        kernels = cell(nNeurons,1);
        min_sp = 3;    % find spikes resulting in transients above min_sp x noise level

        for k = 1:nNeurons
            [C_dec(k,:),S(k,:),kernels{k}] = deconvCa(F_dff(k,:), [], min_sp, true, false, [], 20, [], 0);
        end

        tmp_neuron.S = S;
        tmp_neuron.C_dec = C_dec;
        tmp_neuron.F_dff = F_dff;
        tmp_neuron.F0 = F0;
        
        neurons_f.(['Plane_' num2str(z)]) = tmp_neuron;

    end

    savefast([Filepath '\PospP_CNMF.mat'], 'neurons_f');

end


%%% Old Code

%         % Now find the path to the data, which is needed for refinement
%         [filepath,~,~] = fileparts(cell2mat(Data{i}));
%         One_up = strfind(filepath, '\');
%         
%         % refine spatial components
%         options.spatial_method = 'constrained';
%         [A,b,C,P] = update_spatial_components(Dat, C, f, A, P, options);
% 
%         % refine temporal components
% 
%         [C2,f2,P2,S2,YrA2] = update_temporal_components(Dat,A,b,C,f,P,options);
%         
%         [ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(Dat,A,C,b,f,[],options);
% 
%         if vis
%             % a simple GUI
%             Coor = plot_contours(A,Cn,options,1); close;
%             run_GUI = true;
%             if run_GUI
%                 GUIout = ROI_GUI(M_nr, A, P, options, Cn, C, b, f, Coor, S);   
%                 options = GUIout{2};
%                 keep = GUIout{3};    
%             end
% 
%             % view contour plots of selected and rejected components (optional)
%             throw = ~keep;
%             figure;
%                 ax1 = subplot(121); plot_contours(A(:,keep),Cn,options,0,[],Coor,1,find(keep)); title('Selected components','fontweight','bold','fontsize',14);
%                 ax2 = subplot(122); plot_contours(A(:,throw),Cn,options,0,[],Coor,1,find(throw));title('Rejected components','fontweight','bold','fontsize',14);
%                 linkaxes([ax1,ax2],'xy')
%             % inspect components
% 
%             plot_components_GUI(M_nr,A(:,keep),C(keep,:),b,f,Cn,options);
%         end
%        
%        A_keep = A(:,keep);
%        C_keep = C(keep,:);
%        P2.b = P.b(keep);
%        P2.c1 = P.c1(keep);
%        P2.gn = P.gn(keep);
%        P2.neuron_sn = P.neuron_sn(keep);
%        P2.rval_space = P.rval_space(keep);
%        P2.rval_time = P.rval_time(keep);
%        P2.sn = P.sn;
%        P2.p = P.p;
%        
%        Centroids = Centroids(keep,:);
%         
%        data_lin = reshape(data, [size(data, 2)*size(data, 1), size(data, 3)]);
%        
%         % refine spatial components
%         options.spatial_method = 'constrained';
%         [A,b,C,P] = update_spatial_components(double(data_lin), C_keep, f, A_keep, P2, options);
% 
%         % refine temporal components
% 
%         [C2,f2,P2,S2,YrA2] = update_temporal_components(double(data_lin),A_keep,b,C_keep,f,P2,options);
% 
% 
%         % detrend fluorescence and extract DF/F values
% 
%         options.df_window = 1000; 
%         [F_dff,F0] = detrend_df_f(A,b,C2,f2,YrA2,options);
% 
%         % deconvolve data
% 
%         nNeurons = size(F_dff,1);
%         C_dec = zeros(size(F_dff));
%         S = zeros(size(F_dff));
%         kernels = cell(nNeurons,1);
%         min_sp = 3;    % find spikes resulting in transients above min_sp x noise level
% 
%         for i = 1:nNeurons
%             [C_dec(i,:),S(i,:),kernels{i}] = deconvCa(F_dff(i,:), [], min_sp, true, false, [], 20, [], 0);
%         end
% 
%         % plot a random component
%         if vis
%             figure;imagesc(reshape(max(A_keep, [], 2), [256 256]))
%             figure;imagesc(F_dff, [0 1]); box off;set(gca, 'Ytick', []); xlabel('Time'); ylabel ('# Cells');colormap bone
% 
%             k = randi(nNeurons);
% 
%             figure;plot(F_dff(k,:),'--k'); hold all; plot(C_dec(k,:),'r','linewidth',2);
%                 spt = find(S(k,:));
%                 if spt(1) == 1; spt(1) = []; end
%                 hold on; scatter(spt,repmat(-0.25,1,length(spt)),'m*')
%                 title(['Component ',num2str(k)]);
% 
%                 legend('Fluorescence DF/F','Deconvolved','Spikes')
%         end  
% 
%         Result.S = S;
%         Result.C_dec = C_dec;
%         Result.F_dff = F_dff;
%         Result.F0 = F0;
%         Result.Coor = Coor;  
% 
%         Results{i} = Result;
%         Binary_masks{i} = full(A);
%         Names = [filepath(One_up(end)+1:end) '_' filepath(One_up(end-1)+1:One_up(end)-1)];
% 
% 
%         DataS.(Names).Results = Results;
%         DataS.(Names).Binary_masks = Binary_masks;
% 


