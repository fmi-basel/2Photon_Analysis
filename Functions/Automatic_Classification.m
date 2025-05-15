%% Automatic_Classification
% Code written by Julian Hinz @Luthi lab 2020
% Based on the work of Tianlin Liu

function Automatic_Classification(Filepath)
    
    %Load the trained classifier
    load('OptimalSVM.mat');
    
    % Load dataset
    try
        load([Filepath '\Results_CNMF.mat'])
    catch
       disp('File was not processed, skipping !')
       return
    end
    
    % Next load the Masterstructure containing the Predictors
    try
        load([Filepath '\Master_structure.mat'])
    catch
        disp('MasterStructure not available, creating it !')
        load([Filepath '\Timing_Information.mat'])
        load([Filepath '\Corr_Imgs.mat'])
        Master_structure = Extract_Components(neuron, obj, Timing_Info.FS);
        savefast([Filepath '\Master_structure.mat'], 'Master_structure')
    end
    
    Planes = fieldnames(Master_structure);
        
    % Create new variable for saving classified results
    neurons_f = struct;
    
    for x = 1:size(Planes, 1)
        Cur_Plane = Master_structure.(Planes{x});
        
        % Calculate Distance from Center
        Size = size(Cur_Plane.Corr_imgs);
        Distance_from_center = pdist2([Size(1)/2 Size(2)/2], ...
            Cur_Plane.Centroids) / Size(1);
        
        % Create Predictor Matrix
        PredMatrix = [Cur_Plane.Mean_2_std Cur_Plane.Extracted_values' ...
            Cur_Plane.Area Cur_Plane.circularities Cur_Plane.Perimeter ...
            Cur_Plane.Eccentricity Distance_from_center' ...
            Cur_Plane.Std2crossings Cur_Plane.Correlations' Cur_Plane.AbsoluteSpikeCount ...
            Cur_Plane.Dec_Rawcorr' Cur_Plane.CorrImgs2MedianCorr' Cur_Plane.Sizerel2max ...
            Cur_Plane.Corr_2Local Cur_Plane.Overlaps];
        
        PredMatrix(isinf(PredMatrix)) = 1;
        
        % Run Classifier 
        [Classification, Classification_Score] = predict(OptimalSVM, PredMatrix);
        Delete = logical(Classification);

        % Apply Classifier to the Dataset and create new variable only with
        % neurons that are deemed acceptable
        [tmp_neuron] = update_CNMF_Results(neuron.(['Plane_' num2str(x)]), Delete);
        
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
        tmp_neuron.Classification_Score = Classification_Score;
        
        % Now save all of this in new variable
        neurons_f.(['Plane_' num2str(x)]) = tmp_neuron;
    end
    
    savefast([Filepath '\PospP_CNMF_SVM.mat'], 'neurons_f');
    
end



 