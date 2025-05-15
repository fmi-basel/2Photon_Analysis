%% Get Stats ready from different Datasets to automatically sort components
function Master_structure = Extract_Components(neuron, obj, FS)

% Loop through planes and start extracting statistics about components for
% later sorting

Number_Planes = fieldnames(neuron);
Number_Corr = fieldnames(obj);

    for zz = 1:size(Number_Planes, 1)

        neuron.(Number_Planes{zz}).options.df_window = 1000; 
        Correlation_Imgs = obj.(Number_Corr{zz}).Corr_Img(:,:,4);

        % First measure based on COrr Imgs

        Correlation_Imgs_lin = repmat(reshape(Correlation_Imgs, size(Correlation_Imgs, 1) ...
            * size(Correlation_Imgs, 2), 1), 1, size(neuron.(Number_Planes{zz}).A , 2));

        Extracted_values = sum(Correlation_Imgs_lin .* double(neuron.(Number_Planes{zz}).A > 0), 1)...
            ./sum(neuron.(Number_Planes{zz}).A > 0);

        % Second measure based on Corr Imgs 
        Correlation_Imgs2 = obj.(Number_Corr{zz}).Corr_Img(:,:,1);

        Correlation_Imgs_lin2 = repmat(reshape(Correlation_Imgs2, size(Correlation_Imgs2, 1) ...
            * size(Correlation_Imgs2, 2), 1), 1, size(neuron.(Number_Planes{zz}).A , 2));

        Extracted_values2 = sum(Correlation_Imgs_lin2 .* double(neuron.(Number_Planes{zz}).A > 0), 1)...
            ./sum(neuron.(Number_Planes{zz}).A > 0);

        Relative_Corr = Extracted_values2./median(Correlation_Imgs2(:));
        Norm_Rel_Corr = Relative_Corr./max(Relative_Corr(:));

        % SNR
        Mean_2_std = mean(movmean(zscore(neuron.(Number_Planes{zz}).C), FS*2, 2), 2) ...
            ./std(movmean(zscore(neuron.(Number_Planes{zz}).C), FS*2, 2), [], 2);

        % Spatial Components
        Size = size(neuron.(Number_Planes{zz}).Cn);
        Test_area = reshape(neuron.(Number_Planes{zz}).A, Size(1), Size(2), ...
            size(neuron.(Number_Planes{zz}).A, 2)) > 0;

        Centroids = nan(size(Test_area, 3), 2);
        Area = nan(size(Test_area, 3), 1);
        Perimeter = nan(size(Test_area, 3), 1);
        Eccentricity = nan(size(Test_area, 3), 1);

        parfor z = 1:size(Test_area, 3)
            try
                stats = regionprops(Test_area(:,:,z), 'Centroid', ...
                    'Eccentricity', 'Perimeter', 'PixelIdxList', 'Area');
                Centroids(z, :) = round(stats(1).Centroid);
                Area(z) = stats(1).Area/sqrt((Size(1)*Size(2)));
                Perimeter(z) = stats(1).Perimeter/sqrt((Size(1)*Size(2)));
                Eccentricity(z) = stats(1).Eccentricity;
            end

        end

        distance_from_centre = pdist2([Size(1)/2 Size(2)/2], Centroids) / Size(1);
        circularities = (4 * pi * Area) ./ Perimeter .^2;

        % Now calculate activity metrics
        [F_dff, ~] = detrend_df_f(neuron.(Number_Planes{zz}).A,neuron.(Number_Planes{zz}).b, ...
        neuron.(Number_Planes{zz}).C,neuron.(Number_Planes{zz}).f, ...
        neuron.(Number_Planes{zz}).YrA,neuron.(Number_Planes{zz}).options);

        nNeurons = size(F_dff,1);
        C_dec = zeros(size(F_dff));
        S = zeros(size(F_dff));
        kernels = cell(nNeurons,1);
        min_sp = 3;    % find spikes resulting in transients above min_sp x noise level

        for k = 1:nNeurons
            try
                [C_dec(k,:),S(k,:), kernels{k}] = deconvCa(F_dff(k,:), [], min_sp, true, false, [], 20, [], 0);
                disp(['Finished Neuron ' num2str(k)])
            end
        end

        Correlation_Dec = arrayfun(@(k) corr(C_dec(k,:)', F_dff(k,:)'), 1:size(F_dff, 1), 'Uni', 1);

        % Exclude cells based on their correlation and distance 

        U = logical(triu(ones(size(neuron.(Number_Planes{zz}).Centroids, 1), ...
            size(neuron.(Number_Planes{zz}).Centroids, 1))));
        Distances = zeros(size(neuron.(Number_Planes{zz}).Centroids, 1), ...
            size(neuron.(Number_Planes{zz}).Centroids, 1));

        for kk = 1:size(Centroids, 1)
            Distances(U(kk, :), kk) = round(mean(abs(repmat(Centroids(kk, :), sum(U(kk, :)), 1) ...
                - Centroids(U(kk, :), :)), 2));
        end

        % Derive measure for correlation
        Correlations = corrcoef(movmean(neuron.(Number_Planes{zz}).C, 10, 2)');
        Correlations(find(eye(size(Correlations, 2)))) = 0;
        Max_Correlation = max(Correlations);

        % Size relative to max
        Sizerel2max = Area./max(Area(:));

        % Correlation within ROI relative to local correlation
        Correlation_Imgs2 = obj.(Number_Corr{zz}).Corr_Img(:,:,1);
        LocalArea = round(50*(512/size(Correlation_Imgs2, 1)));
        Corr_2Local = nan(size(Centroids, 1), 1);

        parfor z = 1:size(Centroids, 1)

            % Find Bounding Box
            Area_tmp = [Centroids(z, :) - LocalArea/2; ...
                Centroids(z, :) + LocalArea/2];           
            Area_tmp(Area_tmp < 1) = 1;
            Area_tmp(Area_tmp > size(Correlation_Imgs2, 1)) = ...
                size(Correlation_Imgs2, 1);       

            %Calculate metric
            LocalCorr = Correlation_Imgs2(Area_tmp(1,1):Area_tmp(2,1), ...
                Area_tmp(1,2):Area_tmp(2,2));
            Corr_2Local(z) = Extracted_values2(z)/mean(LocalCorr(:));

        end

        % Overlap with other ROIS

        Num_Neurons = size(neuron.(Number_Planes{zz}).A, 2);
        Spatial_Components = neuron.(Number_Planes{zz}).A;
        Weighted_Components = double(Spatial_Components./max(Spatial_Components) > 0.5);
        Overlaps = nan(Num_Neurons, 1);

        parfor z = 1:Num_Neurons
            Size_Comp = sum(Weighted_Components(:, z)); % Size of component
            Temp_comp = Weighted_Components - repmat(Weighted_Components(:, z), 1, Num_Neurons); % subtract location of the component in question from all the others
            Tmp_Sum = abs(sum(Temp_comp < 0) - Size_Comp);
            Tmp_Sum(z) = 0;
            Overlaps(z) = max(Tmp_Sum / Size_Comp); % Now 'normalize' for size    
        end

        % Number of transients
        Mean2Std = mean(movmean(F_dff, FS*2, 2), 2) + 2*std(movmean(F_dff, FS*2, 2), [], 2);
        Std_Crossings = cell2mat(arrayfun(@(k) gt(movmean(F_dff(k,:), FS*2, 2), Mean2Std(k)), ...
            1:size(F_dff, 1), 'Uni', 0)');
        Number_of_Crossings = nan(Num_Neurons, 1);

        parfor z = 1:Num_Neurons
            [Cross_On, Cross_OFF] = find_sequences_Licks(Std_Crossings(z, :), 1);
            Number_of_Crossings(z) = sum(Cross_OFF - Cross_On > 4);
        end

        % Save all the variables
        Master_structure.(Number_Planes{zz}).Correlations = Max_Correlation;
        Master_structure.(Number_Planes{zz}).Distances = Distances;
        Master_structure.(Number_Planes{zz}).C_dec = C_dec;
        Master_structure.(Number_Planes{zz}).S = S;
        Master_structure.(Number_Planes{zz}).F_dff = F_dff;
        Master_structure.(Number_Planes{zz}).Eccentricity = Eccentricity;
        Master_structure.(Number_Planes{zz}).Perimeter = Perimeter;
        Master_structure.(Number_Planes{zz}).circularities = circularities;
        Master_structure.(Number_Planes{zz}).Area = Area;
        Master_structure.(Number_Planes{zz}).Centroids = Centroids;
        Master_structure.(Number_Planes{zz}).Mean_2_std = Mean_2_std;
        Master_structure.(Number_Planes{zz}).Extracted_values = Extracted_values;
        Master_structure.(Number_Planes{zz}).Corr_imgs = obj.(Number_Corr{zz}).Corr_Img(:,:,4);
        Master_structure.(Number_Planes{zz}).A = neuron.(Number_Planes{zz}).A;
        Master_structure.(Number_Planes{zz}).Bound = neuron.(Number_Planes{zz}).Boundries;
        Master_structure.(Number_Planes{zz}).distance_from_centre = distance_from_centre;

        Master_structure.(Number_Planes{zz}).AbsoluteSpikeCount = sum(S, 2);
        Master_structure.(Number_Planes{zz}).Dec_Rawcorr = Correlation_Dec;
        Master_structure.(Number_Planes{zz}).CorrImgs2MedianCorr = Norm_Rel_Corr;
        Master_structure.(Number_Planes{zz}).Sizerel2max = Sizerel2max;
        Master_structure.(Number_Planes{zz}).Corr_2Local = Corr_2Local;
        Master_structure.(Number_Planes{zz}).Overlaps = Overlaps;
        Master_structure.(Number_Planes{zz}).Std2crossings = Number_of_Crossings;
        
    end

end
