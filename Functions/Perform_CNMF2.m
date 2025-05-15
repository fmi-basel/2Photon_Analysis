%% Implement the Paninksi Workflow for our Data

function Perform_CNMF2(File_Path, Max_Runtime, FS)

    % Check for previous processing reults and prevent overwrite by
    % returning to main Function
    if exist([File_Path '\Results_CNMF.mat']) ~= 0
        disp('File was already processed, skip to the next processing step !')
        return            
    end
        
    % Load data and get sizes         
    Data_handle = matfile([File_Path '\DS_Dat.mat']);     
    Sizes = size(Data_handle.DS_Dat);    
   
    % Set basic parameters of the CNMF
    sizY = Sizes(1:3);
    patch_size = [50,50];                   % size of each patch along each dimension (optional, default: [32,32])
    overlap = [8,8];                        % amount of overlap in each dimension (optional, default: [4,4])

    patches = construct_patches(sizY(1:end-1),patch_size,overlap);

    if sizY(1) < 270       
        tau = 3; % std of gaussian kernel (half size of neuron) 
    elseif sizY(1) < 540     
        tau = 5;  % std of gaussian kernel (half size of neuron) 
    end
    
    p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
    K = [];
    
    if numel(Sizes) <= 3
        Iterations = 1;
    else
        Iterations = Sizes(4);
    end
    
    % Find Optimal minimum Correlation
    obj = [];
    load([File_Path '\Corr_Imgs.mat'])
    
    minCorrelation_value = nan(Iterations, 4);
    for x = 1:Iterations
        Corr_tmp = obj.(['Cor_Images_' num2str(x)]).Corr_Img(:,:,1);
        minCorrelation_value(x, 5) = mean(Corr_tmp(:))+ 3*std(Corr_tmp(:));
        minCorrelation_value(x, 2) = mean(Corr_tmp(:)) + 2*std(Corr_tmp(:));
        minCorrelation_value(x, 3) = mean(Corr_tmp(:))+ 1*std(Corr_tmp(:));
        minCorrelation_value(x, 4) = mean(Corr_tmp(:));
        minCorrelation_value(x, 1) = median(Corr_tmp(:))/2;
    end
    
    completed = zeros(Iterations, 1);

    for k = 1:Iterations
        
        start = tic;
              
        % Load Data
        if Iterations > 1
            data = Data_handle.DS_Dat(:,:,:,k);
        else
            data = Data_handle.DS_Dat;
        end
        
        % run CNMF algorithm on patches and combine - Since there were
        % problems with processing individual planes,try to run the Code as
        % long as the plane is not processed or exit after a while and then
        % I need to check what is wrong with the data
        
        Finished_processing = false;
        counter = 1;
        
        
        while ~Finished_processing
            

            options = CNMFSetParms(...
                'd1',sizY(1),'d2',sizY(2),...
                'temporal_iter',2,...                       % number of block-coordinate descent steps 
                'init_method', 'greedy_corr',...  
                'fr', FS, ...
                'ssub', 1,...                                % downsample in space
                'tsub', round(FS*2),...                 % downsample in time
                'merge_thr', 0.8,...                         % merging threshold
                'min_corr', minCorrelation_value(k, counter),...
                'gnb', 1,...                                 % number of background components
                'spatial_method','regularized',...
                'min_size_thr', floor(10/(540/sizY(1)))); 

            try
                [A, b, C, f, ~, P, ~, YrA] = run_CNMF_patches(data,K,patches,tau,p,options);
                Finished_processing = true;
            catch
                fprintf('\n\n###################################################################\n')
                fprintf(' Error processing Plane - Retrying !!\n')
                fprintf('###################################################################\n\n')
                counter = counter + 1;
            end
            
            if counter > 5
                fprintf('\n\n###################################################################\n')
                fprintf(' Maximum Runtime exceeded - Exiting analysis !!\n')
                fprintf('###################################################################\n\n')
                fprintf('\n')
                error(['CaImAn failed to process plane: ' [File_Path '\DS_Dat'] ...
                    ' Plane: ' num2str(k) ', please check whats up with this file !'])
                fprintf('\n')
            end
        end
        
        Cn = obj.(['Cor_Images_' num2str(k)]).Corr_Img(:,:,1);
        
        [ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep_n] = ...
            classify_components(data, A, C, b, f, [], options);
        
        % Also remove all border cells an
        Outer_area = 5;
        nums = size(A, 2);
        ROI_Shape = reshape(full(A), size(Cn, 1), size(Cn, 2), nums);

        Exclude = ~any([reshape(ROI_Shape(1:Outer_area, : , :), Outer_area*size(Cn, 2), nums); ...
            reshape(ROI_Shape(end-Outer_area+1:end,: , :), Outer_area*size(Cn, 2), nums); ...
            reshape(ROI_Shape(:, 1:Outer_area, :), Outer_area*size(Cn, 1), nums); ...
            reshape(ROI_Shape(:, end-Outer_area+1:end, :), Outer_area*size(Cn, 1), nums)] > 0);
        keep = Exclude'& keep_n;
        
        % Apply the changes that come with the classification
        
        A = A(:,keep);
        C = C(keep,:);
        YrA =  YrA(keep,:);
        
        P.b = P.b(keep);
        P.c1 = P.c1(keep);
        P.gn = P.gn(keep);
        P.neuron_sn = P.neuron_sn(keep);
        P.rval_space = P.rval_space(keep);
        P.rval_time = P.rval_time(keep);
        P.sn = P.sn;
        P.p = P.p;

        try
            Coor = plot_contours(A, Cn, options, 1);close;
        catch
            disp('Hey')
        end
        
        Centroids = nan(size(Coor, 1), 2);
        Delete = zeros(size(Coor, 1), 1);
        Boundries = cell(size(Coor, 1), 1);
        
        for kk = 1:size(Coor, 1)
            Bin_Image = zeros(sqrt(size(A, 1)), sqrt(size(A, 1)));
            Indices =  transpose(Coor{kk});
            if isempty(Indices)
                Delete(kk) = 1;
                continue
            end
            Bin_Image(sub2ind([sqrt(size(A, 1)), sqrt(size(A, 1))], Indices(:,1), ...
                Indices(:,2))) = 1;
            stats = regionprops(Bin_Image, 'Centroid');
            Boundries{kk} = bwboundaries(Bin_Image); 
            
            Centroids(kk, :) = round(stats.Centroid);
        end
        
        A = full(A);
        
        % Exclude cells based on their correlation and distance 
        
        U = logical(triu(ones(size(Centroids, 1), size(Centroids, 1))));
        Distances = zeros(size(Centroids, 1), size(Centroids, 1));

        for kk = 1:size(Centroids, 1)
            Distances(U(kk, :), kk) = round(mean(abs(repmat(Centroids(kk, :), sum(U(kk, :)), 1) ...
                - Centroids(U(kk, :), :)), 2));
        end

        Correlations = corrcoef(movmean(C, 10, 2)');

        % Merge Neurons

        Merge_thr(1) = 0.6; Merge_thr(2) = 10; 
        
        Corr_Dist_del = Correlations > Merge_thr(1) & Distances < Merge_thr(2); %| Distances < 3;
        Corr_Dist_del(U ~= 0) = 0;

        [Row, Column] = ind2sub(size(Correlations), find(Corr_Dist_del == 1));
        Merge_DelInd_Sub = [Row Column];

        % Collect all the results in a structure
        
        neuron.(['Plane_' num2str(k)]).A = A;
        neuron.(['Plane_' num2str(k)]).b = b;
        neuron.(['Plane_' num2str(k)]).C = C;
        neuron.(['Plane_' num2str(k)]).f = f;
        neuron.(['Plane_' num2str(k)]).P = P;
        neuron.(['Plane_' num2str(k)]).options = options;
        neuron.(['Plane_' num2str(k)]).Delete = Delete;
        neuron.(['Plane_' num2str(k)]).patches = patches;
        neuron.(['Plane_' num2str(k)]).Centroids = Centroids;
        neuron.(['Plane_' num2str(k)]).Cn = Cn;
        neuron.(['Plane_' num2str(k)]).Coor = Coor;
        neuron.(['Plane_' num2str(k)]).Merge_DelInd_Sub = Merge_DelInd_Sub;
        neuron.(['Plane_' num2str(k)]).Distances = Distances;
        neuron.(['Plane_' num2str(k)]).Boundries = Boundries;
        neuron.(['Plane_' num2str(k)]).YrA = YrA;
        
        One_up = strfind(File_Path, '\');    
        stop = toc(start);
        
        disp('******************************************************')
        fprintf(['Finished with Imaging file ' num2str(k) '/' num2str(Iterations)...
            ' in ' File_Path(One_up(end-2)+1:One_up(end-1)-1) ...
            ' in folder ' File_Path(One_up(end-1)+1:One_up(end)-1) '\n'])
        fprintf(['CNMF took ' num2str(stop) ' s to complete \n'])
        disp('******************************************************')
        
        completed(k) = true;
        
    end
    
    if any(completed)
        % Save the results to process later
        savefast([File_Path '\Results_CNMF.mat'], 'neuron');
        
        % Extract stats from components for later analysis
        Master_structure = Extract_Components(neuron, obj, FS);
        savefast([File_Path '\Master_structure.mat'], 'Master_structure')
        
    end
        
end
    

% 
%         for k = 1:size(Merge_DelInd_Sub, 1)
%             figure;
%             subplot(1,2,1);hold on
%             plot(C(Merge_DelInd_Sub(k,1), :), 'k')
%             plot(C(Merge_DelInd_Sub(k,2), :), 'r')
%             subplot(1,2,2)
%             imagesc(reshape(max(A(:, Merge_DelInd_Sub(k,:)), [], 2), size(Cn, 1), size(Cn, 2)))
%         end

%         manual_select = false;

%         if manual_select
%             figure;
%             scatter(abs(distance_Correlation(:, 1)), distance_Correlation(:, 2), 10, 'filled', 'ok');
%             box off; xlabel('Correlation'); ylabel('Distance [pix]')
%             [Merge_thr(2), Merge_thr(1)] = ginput(1);
%             hold on;
%             rectangle('Position',[Merge_thr(2) 0 1-Merge_thr(2) Merge_thr(1)], 'EdgeColor',  'r', 'linewidth', 2)
%             drawnow
%         end


%         Correlations_lin = reshape(Correlations, [1, size(Correlations, 1)*size(Correlations, 2)]);
%         Correlations_lin_n = Correlations_lin;
%         Correlations_lin_n(U == 1) = 0;
%         distance_Correlation = [transpose(Correlations_lin_n(U ~= 1)) Distances(U ~= 1)];

%         if vis
%             figure;imagesc(Correlations)
% 
%             figure; 
%             subplot(1,2,1); 
%             imagesc(Distances); set(gca, 'Xtick', [], 'Ytick', []); box off; 
%             subplot(1,2,2)
%             imagesc(reshape(Correlations_lin_n, [size(Correlations, 1) size(Correlations, 2)]))
%             set(gca, 'Xtick', [], 'Ytick', []); box off; colormap bone
% 
%             figure;scatter(abs(distance_Correlation(:, 1)), distance_Correlation(:, 2)/size(data, 1), 10, 'filled', 'ok');
%             box off; xlabel('Correlation'); ylabel('Rel. Distance')
%         end
