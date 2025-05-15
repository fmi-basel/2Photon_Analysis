%% Derive stats from manually sorted neurons to derive an automatic classifier 

clear 
close
clc

addpath(genpath('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Github\gluthi-2Photon'));

Directory = '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Imaging_2p\US\Experiments_August2019\Masterstructures\';
Filenames = dir([Directory '*.mat']);
Filenames_Sorted = sort_nat({Filenames.name});

% Struct_all = struct;

Collection_A = cell(size(Filenames, 1), 1);
Collection_F_dff = cell(size(Filenames, 1), 1);
Collection_C_dec = cell(size(Filenames, 1), 1);

Collection_Mean_2_std = cell(size(Filenames, 1), 1);
Collection_Extracted_values = cell(size(Filenames, 1), 1);
Collection_Area = cell(size(Filenames, 1), 1);
Collection_circularities = cell(size(Filenames, 1), 1);
Collection_Perimeter = cell(size(Filenames, 1), 1);
Collection_Eccentricity = cell(size(Filenames, 1), 1);
Collection_Distance_from_center = cell(size(Filenames, 1), 1);

Collection_Std2crossings = cell(size(Filenames, 1), 1);
Collection_Correlations = cell(size(Filenames, 1), 1);
Collection_SpikeCount = cell(size(Filenames, 1), 1);
Collection_Dec_Rawcorrr = cell(size(Filenames, 1), 1);
Collection_CorrImgs2MedianCorr = cell(size(Filenames, 1), 1);
Collection_Sizerel2max = cell(size(Filenames, 1), 1);
Collection_Corr_2Local = cell(size(Filenames, 1), 1);
Collection_Overlaps = cell(size(Filenames, 1), 1);


Collection_Decision = cell(size(Filenames, 1), 1);

for zz = 1:size(Filenames, 1)
    % Load the extracted data
    load([Directory Filenames_Sorted{zz}])
    % Struct_all.(['Ses_' num2str(zz)]) = Master_structure.(['Ses_' num2str(zz)]);
    
    Planes_names = fieldnames(Master_structure.(['Ses_' num2str(zz)]));
    Temp_cell = cell(size(Planes_names, 1) - 1, 19);
    
    for z = 1:size(Planes_names, 1) - 1
        Size = size(Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).Corr_imgs);
        Temp_cell{z, 1} = pdist2([Size(1)/2 Size(2)/2], ...
            Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).Centroids) / Size(1);
        Temp_cell{z, 2} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).Extracted_values;
        Temp_cell{z, 3} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).Area;
        Temp_cell{z, 4} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).circularities;
        Temp_cell{z, 5} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).Perimeter;
        Temp_cell{z, 6} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).Eccentricity;
        Temp_cell{z, 7} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).Mean_2_std;
        Temp_cell{z, 8} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).F_dff;
        Temp_cell{z, 9} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).A;
        Temp_cell{z, 10} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).C_dec;
        Temp_cell{z, 11} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).Decision{2};
        Temp_cell{z, 12} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).Correlations;
        Temp_cell{z, 13} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).AbsoluteSpikeCount;
        Temp_cell{z, 14} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).Dec_Rawcorr;
        Temp_cell{z, 15} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).CorrImgs2MedianCorr;
        Temp_cell{z, 16} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).Sizerel2max;
        Temp_cell{z, 17} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).Corr_2Local;
        Temp_cell{z, 18} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).Overlaps;
        Temp_cell{z, 19} = Master_structure.(['Ses_' num2str(zz)]).(Planes_names{z}).Std2crossings;
        
    end
    
    Collection_A{zz} = cell2mat(Temp_cell(:, 9)');
    Collection_F_dff{zz} = cell2mat(Temp_cell(:, 8));
    Collection_C_dec{zz} = cell2mat(Temp_cell(:, 10));
    
    Collection_Mean_2_std{zz} = cell2mat(Temp_cell(:, 7));
    Collection_Extracted_values{zz} = cell2mat(Temp_cell(:, 2)')';
    Collection_Area{zz} = cell2mat(Temp_cell(:, 3));
    Collection_circularities{zz} = cell2mat(Temp_cell(:, 4));
    Collection_Perimeter{zz} = cell2mat(Temp_cell(:, 5));
    Collection_Eccentricity{zz} = cell2mat(Temp_cell(:, 6));
    Collection_Distance_from_center{zz} = cell2mat(Temp_cell(:, 1)')';
    Collection_Decision{zz} =  cell2mat(Temp_cell(:, 11));
    
    Collection_Std2crossings{zz} = cell2mat(Temp_cell(:, 19));
    Collection_Correlations{zz} = cell2mat(Temp_cell(:, 12)');
    Collection_SpikeCount{zz} = cell2mat(Temp_cell(:, 13));
    Collection_Dec_Rawcorrr{zz} = cell2mat(Temp_cell(:, 14)');
    Collection_CorrImgs2MedianCorr{zz} = cell2mat(Temp_cell(:, 15)');
    Collection_Sizerel2max{zz} = cell2mat(Temp_cell(:, 16));
    Collection_Corr_2Local{zz} = cell2mat(Temp_cell(:, 17));
    Collection_Overlaps{zz} = cell2mat(Temp_cell(:, 18));


    disp(['Finished loading File ' num2str(zz) ' out of ' num2str(size(Filenames, 1))])
end

Mean_2_std = cell2mat(Collection_Mean_2_std);
Extracted_values = cell2mat(Collection_Extracted_values);
Area = cell2mat(Collection_Area);
circularities = cell2mat(Collection_circularities);
Perimeter = cell2mat(Collection_Perimeter);
Eccentricity = cell2mat(Collection_Eccentricity);
Distance_from_center = cell2mat(Collection_Distance_from_center);
Decision = cell2mat(Collection_Decision);

Std2crossings = cell2mat(Collection_Std2crossings);
Correlations = cell2mat(Collection_Correlations')';
SpikeCount = cell2mat(Collection_SpikeCount);
Dec_Rawcorrr = cell2mat(Collection_Dec_Rawcorrr')';
CorrImgs2MedianCorr = cell2mat(Collection_CorrImgs2MedianCorr')';
Sizerel2max = cell2mat(Collection_Sizerel2max);
Corr_2Local = cell2mat(Collection_Corr_2Local);
Overlaps = cell2mat(Collection_Overlaps);


clear Collection_Mean_2_std Collection_Extracted_values Collection_Area Collection_circularities
clear Collection_Perimeter Collection_Eccentricity Collection_Distance_from_center
clear Collection_Std2crossings Collection_Correlations Collection_SpikeCount Collection_Dec_Rawcorrr
clear Collection_CorrImgs2MedianCorr Collection_Sizerel2max Collection_Corr_2Local Collection_Overlaps


Matrix_Cellspace = [Mean_2_std Extracted_values Area circularities Perimeter Eccentricity Distance_from_center ...
    Std2crossings Correlations SpikeCount Dec_Rawcorrr CorrImgs2MedianCorr Sizerel2max Corr_2Local...
    Overlaps];

Matrix_Cellspace(isnan(Matrix_Cellspace)) = 0;
Matrix_Cellspace(isinf(Matrix_Cellspace)) = 1;


Matrix_Cellspace_norm = Matrix_Cellspace./max(Matrix_Cellspace);
[coeff,score,latent,tsquared,explained,mu] = pca(Matrix_Cellspace_norm);

figure; plot(cumsum(explained), 'ok', 'linewidth', 2);
xlim([0 size(Matrix_Cellspace, 2)]); ylim([0 100]);
xlabel('# Components'); ylabel('Variance explained');box off

figure;hold on;
scatter3(score(Decision, 1), score(Decision, 2), score(Decision, 3), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter3(score(~Decision, 1), score(~Decision, 2), score(~Decision, 3), 3, 'filled', 'r')

figure;
scatter3(Matrix_Cellspace_norm(:, 1), Matrix_Cellspace_norm(:, 2), Matrix_Cellspace_norm(:, 7), 3, 'filled', 'k')

figure;hold on;
scatter3(Matrix_Cellspace_norm(Decision, 14), Matrix_Cellspace_norm(Decision, 13), Matrix_Cellspace_norm(Decision, 11), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter3(Matrix_Cellspace_norm(~Decision, 14), Matrix_Cellspace_norm(~Decision, 13), Matrix_Cellspace_norm(~Decision, 11), 3, 'filled', 'r')


figure;
hold on;
[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1), Matrix_Cellspace_norm(:, 1), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 1), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 1), 3, 'filled', 'k')

[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 1, Matrix_Cellspace_norm(:, 2), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 2), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 2), 3, 'filled', 'b')

[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 2, Matrix_Cellspace_norm(:, 3), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 3), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 3), 3, 'filled', 'r')

[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 3, Matrix_Cellspace_norm(:, 4), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 4), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 4), 3, 'filled', 'g')

[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 4, Matrix_Cellspace_norm(:, 5), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 5), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 5), 3, 'filled', 'm')

[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 5, Matrix_Cellspace_norm(:, 6), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 6), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 6), 3, 'filled', 'c')

[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 6, Matrix_Cellspace_norm(:, 7), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 7), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 7), 3, 'filled', 'y')

[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 7, Matrix_Cellspace_norm(:, 8), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 8), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 8), 3, 'filled', 'k')

[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 8, Matrix_Cellspace_norm(:, 9), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 9), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 9), 3, 'filled', 'b')

[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 9, Matrix_Cellspace_norm(:, 10), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 10), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 10), 3, 'filled', 'r')

[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 10, Matrix_Cellspace_norm(:, 11), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 11), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 11), 3, 'filled', 'g')

[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 11, Matrix_Cellspace_norm(:, 12), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 12), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 12), 3, 'filled', 'm')

[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 12, Matrix_Cellspace_norm(:, 13), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 13), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 13), 3, 'filled', 'c')

[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 13, Matrix_Cellspace_norm(:, 14), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 14), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 14), 3, 'filled', 'y')

[X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 14, Matrix_Cellspace_norm(:, 15), 0.5);
scatter(X(Decision), Matrix_Cellspace_norm(Decision, 15), 3, 'filled', 'MarkerFaceColor', [0.7 0.7,0.7])
scatter(X(~Decision), Matrix_Cellspace_norm(~Decision, 15), 3, 'filled', 'r')

xlim([0 16]); box off;ylim([0 1])
set(gca, 'XTick', [1:15], 'XTickLabel', {'SNR', 'Signal Strength', 'Area', ...
    'circularity', 'Perimeter', 'Eccentricity', 'Distance from Center', 'STDCrossings' ...
    'MaxCorrelation', 'SumSpikes', 'CorrDecRaw', 'LocCor2MedTotalCorr', 'Size2max', ...
    'LocalCorr2ROI', 'Overlap'}, 'FontSize', 24);
ylabel('norm. variables', 'FontSize', 24)
xtickangle(45)

% 
% figure;
% hold on;
% [X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1), Matrix_Cellspace_norm(:, 1), 0.5);
% scatter(X, Matrix_Cellspace_norm(:, 1), 3, 'filled', 'k')
% 
% [X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 1, Matrix_Cellspace_norm(:, 2), 0.5);
% scatter(X, Matrix_Cellspace_norm(:, 2), 3, 'filled', 'b')
% 
% [X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 2, Matrix_Cellspace_norm(:, 3), 0.5);
% scatter(X, Matrix_Cellspace_norm(:, 3), 3, 'filled', 'r')
% 
% [X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 3, Matrix_Cellspace_norm(:, 4), 0.5);
% scatter(X, Matrix_Cellspace_norm(:, 4), 3, 'filled', 'g')
% 
% [X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 4, Matrix_Cellspace_norm(:, 5), 0.5);
% scatter(X, Matrix_Cellspace_norm(:, 5), 3, 'filled', 'm')
% 
% [X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 5, Matrix_Cellspace_norm(:, 6), 0.5);
% scatter(X, Matrix_Cellspace_norm(:, 6), 3, 'filled', 'c')
% 
% [X] = Viola_Plot_taken(ones(size(Matrix_Cellspace, 1), 1) + 6, Matrix_Cellspace_norm(:, 7), 0.5);
% scatter(X, Matrix_Cellspace_norm(:, 7), 3, 'filled', 'y')
% 
% xlim([0 8]); box off
% set(gca, 'XTick', [1:7], 'XTickLabel', {'SNR', 'Signal Strength', 'Area', ...
%     'circularity', 'Perimeter', 'Eccentricity', 'Distance from Center'});
% ylabel('norm. variables')
% xtickangle(45)




%%% SVM classifier

Dataset = Matrix_Cellspace_norm;
Class_label = Decision;

Mdl = fitcsvm(Dataset,Class_label,'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'));


savefast('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Imaging_2p\US\Experiments_August2019\ML\Classlabels.mat', 'Class_label')
savefast('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Imaging_2p\US\Experiments_August2019\ML\Dataset.mat', 'Dataset')




%%%%% 
% Old Code
% Correlation_Imgs = obj.Cor_Images_1.Corr_Img(:,:,4);
% ROIS = neuron.Plane_1.A;
% Correlation_Imgs_lin = repmat(reshape(Correlation_Imgs, size(Correlation_Imgs, 1) ...
%     * size(Correlation_Imgs, 2), 1), 1, size(ROIS, 2));
% 
% Extracted_values = sum(Correlation_Imgs_lin .* double(ROIS > 0), 1)./sum(ROIS > 0);
% figure; hist(Extracted_values./sum(ROIS > 0), [0:0.05:5])
% 
% 
% 
% 
% g = figure;
% imagesc(Correlation_Imgs); 
% colormap gray;hold on;
% for mm = 1:size(neuron.Plane_1.Boundries, 1)
%     cont = neuron.Plane_1.Boundries{mm}; 
%     if Neurons_tocheck(mm)
%         plot(cont{1}(:,1),cont{1}(:,2),'Color', 'b', 'linewidth', 0.75); 
%     else
%         plot(cont{1}(:,1),cont{1}(:,2),'Color', 'r', 'linewidth', 0.75); 
%     end
% end
% 
% figure; hist(mean(neuron.Plane_1.C, 2)./std(mean(neuron.Plane_1.C, 2)));
% 
% Selection = mean(neuron.Plane_1.C, 2)./std(mean(neuron.Plane_1.C, 2)) > 0.6;
% Neurons_tocheck = Extracted_values./sum(ROIS > 0) > 0.4;
% Size_crit = sum(ROIS > 0) > 50;
% 
% g = figure;
% imagesc(Correlation_Imgs); 
% colormap gray;hold on;
% for mm = 1:size(neuron.Plane_1.Boundries, 1)
%     cont = neuron.Plane_1.Boundries{mm}; 
%     if Selection(mm) && Neurons_tocheck(mm) && Size_crit(mm)
%         plot(cont{1}(:,1),cont{1}(:,2),'Color', 'b', 'linewidth', 0.75); 
%     else
%         plot(cont{1}(:,1),cont{1}(:,2),'Color', 'r', 'linewidth', 0.75); 
%     end
% end
% 
% ROIS_to_plot = find(~Selection & ~Neurons_tocheck');
% figure;
% subplot(1,2,1)
% hold on;
% for z = 1:numel(ROIS_to_plot)
%     plot(zscore(neuron.Plane_1.C(ROIS_to_plot(z), :)) + z*10, 'r')
% end
% ROIS_to_plot = find(Selection & Neurons_tocheck');
% 
% subplot(1,2,2)
% hold on;
% for z = 1:numel(ROIS_to_plot)
%     plot(zscore(neuron.Plane_1.C(ROIS_to_plot(z), :)) + z*10, 'b')
% end
% 
% 
% 
% 
% 
% ROIS_to_plot = find(~Selection & ~Neurons_tocheck');
% for z = 1:numel(ROIS_to_plot)
%     figure
%     plot(zscore(neuron.Plane_1.C(ROIS_to_plot(z), :)) + z*10, 'b')
%     title(['Neuron ' num2str(ROIS_to_plot(z))])
% end
% 
% 
% 
% ROIS2 = neurons_f.Plane_1.A;
% Correlation_Imgs_lin = repmat(reshape(Correlation_Imgs, size(Correlation_Imgs, 1) ...
%     * size(Correlation_Imgs, 2), 1), 1, size(ROIS2, 2));
% 
% Extracted_values2 = max(Correlation_Imgs_lin .* double((ROIS2 > 0)),[], 1);
% figure; hist(Extracted_values2, [0:0.05:5])
% 
% Neurons_tocheck = Extracted_values2 < 1;
% 
% g = figure;
% imagesc(Correlation_Imgs); 
% colormap gray;hold on;
% for mm = 1:size(neurons_f.Plane_1.Boundries, 1)
%     cont = neurons_f.Plane_1.Boundries{mm}; 
%     if Neurons_tocheck(mm)
%         plot(cont{1}(:,1),cont{1}(:,2),'Color', 'b', 'linewidth', 0.75); 
%     end
% end

