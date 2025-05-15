%% Classification evaluation
% Written by Julian Hinz @Luthi lab 2021

addpath(genpath('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Github\gluthi-2Photon'))

close all
clear all
clc

% Training set definition
TrainingFraction = 0.8;

% Load the Matrix containing the predictors
load('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Classification\TwoPhoton\ClassMatrix.mat')

% Now select part of the data that has a label
Predictors = CompleteMatrix(~isnan(CompleteMatrix(:, 17)), :);
Predictors(isinf(Predictors)) = 1;

% Now split in test and train dataset
Inds_Train = datasample(1:size(Predictors, 1), round(numel(1:size(Predictors, 1))...
    *TrainingFraction), 'Replace', false);
Inds_Test = setdiff(1:size(Predictors, 1), Inds_Train);

% Create the datasets
X_train_raw = Predictors(Inds_Train, 2:16);
X_test_raw = Predictors(Inds_Test, 2:16);

% Seperate indices for later comparison
X_train_raw_Indices = Predictors(Inds_Train, 1);
X_test_raw_Indices = Predictors(Inds_Test, 1);

% Label 0 indicates real neuron; 
y_train_raw = Predictors(Inds_Train, end);
y_test_raw = Predictors(Inds_Test, end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Now Train the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Increments = 0.5;
precision_train_coll = nan(numel(0.5:Increments:10), 9, 2);
count = 1;

SaveModel = cell(numel(0.5:Increments:10), 1);

% SVM
for x = 0.5:Increments:10
    
    SVMModel = fitcsvm(X_train_raw,y_train_raw, 'Standardize', true, 'KernelFunction',...
        'RBF','KernelScale','auto', 'Cost', [0, 1; x, 0]);

    [train_output, train_output_score] = predict(SVMModel, X_train_raw);
    precision_train = Evaluate_precision(y_train_raw, train_output);
    precision_train_coll(count, 3:end, 1) = precision_train;
    
    [test_output,test_output_score] = predict(SVMModel, X_test_raw);
    precision_test = Evaluate_precision(y_test_raw, test_output);
    abs_test_score = abs(test_output_score(:, 1));
    
    FP = sum(test_output == 0 & y_test_raw == 1)/size(test_output, 1);
    FN = sum(test_output == 1 & y_test_raw == 0)/size(test_output, 1);
    
    precision_train_coll(count, 1:2, 2) = [FP; FN];
    precision_train_coll(count, 3:end, 2) = precision_test;
    
    SaveModel{count} = SVMModel;
    
    disp(['Round ' num2str(count) ' out of 20 finished!'])
    
    count = count + 1;
end

% Plot accuracy measures
figure; plot([0.5:Increments:10], precision_train_coll(:, 1, 2), 'or'); 
hold on; plot([0.5:Increments:10], precision_train_coll(:, 2, 2), 'ob'); 
legend('FP', 'FN'); xlabel('Bias'); ylabel('Fraction FP and FN'); box off
title('SVM classifcation accuracy')

%% Decide on Bias - Currently it looks like 2.5 is best
Bias = 2;
OptimalSVM = SaveModel{Bias};

% Now use the model to predict the Output again
[test_output,test_output_score] = predict(OptimalSVM, X_test_raw);
precision_test = Evaluate_precision(y_test_raw, test_output);
abs_test_score = abs(test_output_score(:, 1));

% Calculate FP and FN rate
FP = sum(test_output == 0 & y_test_raw == 1)/size(test_output, 1);
FN = sum(test_output == 1 & y_test_raw == 0)/size(test_output, 1);

% Find indices of the missclasified neurons and extract for further
% analysis
Indices_FP = X_test_raw_Indices(test_output == 0 & y_test_raw == 1);
Indices_FN = X_test_raw_Indices(test_output == 1 & y_test_raw == 0);

% Load the neuron data, select the relevant neurons that were wrongly
% assigned and then delete the huge file again

Neurons_Loaded = load(['\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\'...
    'hinzjuli\Data\Classification\TwoPhoton\Neurons.mat']);

Neurons_FP = Neurons_Loaded.CompleteNeurons(Indices_FP, :);
Neurons_FN = Neurons_Loaded.CompleteNeurons(Indices_FN, :);

clear Neurons_Loaded

%%% Now start working with the wrongly assigned neurons
% First, lets visualize the spatial footprints and Centroids to see if a
% pattern emerges

% Lets plot centroids
Centroids_FP = cell2mat(Neurons_FP(:, 4));
Centroids_FN = cell2mat(Neurons_FN(:, 4));

figure; 
subplot(1,3,1)
scatter(Centroids_FP(:, 1), Centroids_FP(:, 2), 10, 'filled', 'r'); hold on;
scatter(Centroids_FN(:, 1), Centroids_FN(:, 2), 10, 'filled', 'b'); 
box off; xlim([0 256]); ylim([0 256]); set(gca, 'XTick', [], 'Ytick', [])
subplot(1,3,2)
scatter(Centroids_FP(:, 1), Centroids_FP(:, 2), 10, 'filled', 'r'); 
box off; xlim([0 256]); ylim([0 256]); set(gca, 'XTick', [], 'Ytick', [])
subplot(1,3,3)
scatter(Centroids_FN(:, 1), Centroids_FN(:, 2), 10, 'filled', 'b'); 
box off; xlim([0 256]); ylim([0 256]); set(gca, 'XTick', [], 'Ytick', [])

% Lets plot spatial components
SpatialComp_FP = Neurons_FP(:, 1);
SpatialComp_FP = cat(3, SpatialComp_FP{:});
SpatialComp_FN = Neurons_FN(:, 1);
SpatialComp_FN = cat(3, SpatialComp_FN{:});

figure; 
subplot(1,2,1)
imagesc(max(SpatialComp_FP, [], 3)); colormap gray
box off; xlim([0 256]); ylim([0 256]); set(gca, 'XTick', [], 'Ytick', [])
subplot(1,2,2)
imagesc(max(SpatialComp_FN, [], 3)); colormap gray
box off; xlim([0 256]); ylim([0 256]); set(gca, 'XTick', [], 'Ytick', [])

% Next, let's look at the temporal component and see if there is a
% recognizable pattern there

%%% FP - Raw
TemporalCompRaw_FP = Neurons_FP(:, 5);
FS_FP = Neurons_FP(:, 6);

figure; hold on;
for x = 1:size(TemporalCompRaw_FP, 1)
    plot([1:size(TemporalCompRaw_FP{x}, 2)]./FS_FP{x}, ...
        zscore(TemporalCompRaw_FP{x}) + (x-1)*5, 'linewidth', 1, 'color', ...
        rand(1, 3));
end

set(gca, 'Ytick', [])
xlabel('Time [s]'); ylabel('# Neurons')


%%% FP - Decovolved
TemporalCompRaw_FP = Neurons_FP(:, 2);
FS_FP = Neurons_FP(:, 6);

figure; hold on;
for x = 1:size(TemporalCompRaw_FP, 1)
    plot([1:size(TemporalCompRaw_FP{x}, 2)]./FS_FP{x}, ...
        zscore(TemporalCompRaw_FP{x}) + (x-1)*5, 'linewidth', 1, 'color', ...
        rand(1, 3));
end

set(gca, 'Ytick', [])
xlabel('Time [s]'); ylabel('# Neurons')

%%% FN - Raw
TemporalCompRaw_FN = Neurons_FN(:, 5);
[~, indFN] = sort(cellfun(@(c) max(size(c)), TemporalCompRaw_FN), 'descend');
FS_FN = Neurons_FN(:, 6);

figure; hold on;
for x = 1:size(TemporalCompRaw_FN, 1)
    plot([1:size(TemporalCompRaw_FN{indFN(x)}, 2)]./FS_FN{indFN(x)}, ...
        zscore(TemporalCompRaw_FN{indFN(x)}) + (x-1)*5, 'linewidth', 1, 'color', ...
        rand(1, 3));
end

set(gca, 'Ytick', [])
xlabel('Time [s]'); ylabel('# Neurons')
ylim([0 size(TemporalCompRaw_FN, 1)*5])

%%% FN - Deconvolved
TemporalCompRaw_FN = Neurons_FN(:, 2);
[~, indFN] = sort(cellfun(@(c) max(size(c)), TemporalCompRaw_FN), 'descend');
FS_FN = Neurons_FN(:, 6);

figure; hold on;
for x = 1:size(TemporalCompRaw_FN, 1)
    plot([1:size(TemporalCompRaw_FN{indFN(x)}, 2)]./FS_FN{indFN(x)}, ...
        zscore(TemporalCompRaw_FN{indFN(x)}) + (x-1)*5, 'linewidth', 1, 'color', ...
        rand(1, 3));
end

set(gca, 'Ytick', [])
xlabel('Time [s]'); ylabel('# Neurons')
ylim([0 size(TemporalCompRaw_FN, 1)*5])


%%%% Check, why they were excluded, by plotting the average of the
%%%% individual predictors against the average in FN and FP groups

FP_Pred = X_test_raw(test_output == 0 & y_test_raw == 1, :);
FN_Pred = X_test_raw(test_output == 1 & y_test_raw == 0, :);

figure; 
hold on;
for x = 1:size(FP_Pred, 2)
    
    curr = rand(1, 3);
    scatter(zeros(size(FP_Pred, 1), 1) + (x*4) - 2, abs(FP_Pred(:, x)), ...
        20, 'filled','MarkerFaceColor', curr, 'jitter', 'on', ...
        'jitterAmount', 0.3)
     xAxis =x*4-0.5 - 2:0.1:x*4+0.5-2;
     plot(xAxis,  repelem(mean(abs(FP_Pred(:, x)), 1), ...
        numel(xAxis)), 'color', 'k', ...
        'linewidth', 3)
    
    scatter(zeros(size(X_train_raw, 1), 1) + x*4, abs(X_train_raw(:, x)), ...
        20, 'filled','MarkerFaceColor', curr, 'jitter', 'on', ...
        'jitterAmount', 0.3)
     xAxis =x*4-0.5:0.1:x*4+0.5;
     plot(xAxis,  repelem(mean(abs(X_train_raw(:, x)), 1), ...
        numel(xAxis)), 'color', 'k', ...
        'linewidth', 3)

end

set(gca, 'YScale', 'log')
