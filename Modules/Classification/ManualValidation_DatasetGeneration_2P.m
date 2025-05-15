%% Generate a Dataset for manual validation

addpath(genpath('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Github\gluthi-2Photon\'))

clear
close
clc

% First load the data
load('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Classification\TwoPhoton\Neurons.mat')
load('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Classification\TwoPhoton\ClassMatrix.mat')
load('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Github\gluthi-2Photon\Modules\Classification\OptimalSVM.mat')

% Number of components to sort 
Num_Sorting = 3000;

% Generate Predictions
[output,test_output_score] = predict(OptimalSVM, CompleteMatrix(:, 2:17-1));

% Next subselect data with equal sampling in respect to distance from
% decision boundry

Range = [-2.8:0.2:2.8];

figure; 
subplot(1,2,1)
hist(test_output_score(:, 1), Range, 'k'); box off
xlabel('Distance from Decision Boundry')
ylabel('Num Components')
title('All Distances')
xlim([Range(1) Range(end)])
subplot(1,2,2)
idx = ~isnan(CompleteMatrix(:, 17));
hist(test_output_score(idx, 1), Range, 'k'); box off
xlabel('Distance from Decision Boundry')
title('Distances with manual sorting Results')
xlim([Range(1) Range(end)])

% Preallocate - calculate fraction of components that are sorted acording
% to manual label assigned by Julian
Fraction_Match = nan(size(Range, 2) - 1, 1);
test_output_score_sorted = output(idx, 1);
CompleteMatrix_sorted = CompleteMatrix(idx, 17);
test_score_sorted = test_output_score(idx, 1);

for x = 1:size(Range, 2) - 1   
    IndexofRange = find(test_score_sorted(:, 1) > Range(x) & ...
        test_score_sorted(:, 1) < Range(x + 1));   
    Match = CompleteMatrix_sorted(IndexofRange) == 1 & test_output_score_sorted(IndexofRange) == 1 | ...
        CompleteMatrix_sorted(IndexofRange) == 0 & test_output_score_sorted(IndexofRange) == 0;
    Fraction_Match(x) = sum(Match)/numel(IndexofRange)*100;
end

figure; 
bar(Fraction_Match, 'k')
ylabel('Percent sorted corectly according to Julians labels')
xlabel('Distance from Decision Boundry')
box off
StringConcat = [cellfun(@num2str, mat2cell(Range(1:end-1), 1, ones(numel...
    (Range)-1, 1))', 'UniformOutput', false) repelem({'-'}, ...
    numel(Range)-1)' cellfun(@num2str, mat2cell(Range(2:end), 1, ...
    ones(numel(Range)-1, 1))', 'UniformOutput', false)];
StringConcatFull = arrayfun(@(row)[(StringConcat{row,:})],(1:size(...
    StringConcat,1))','UniformOutput', false);
set(gca, 'XTick', 1:size(Range, 2) - 1, 'XTickLabel', StringConcatFull)
xtickangle(45)

% Now generate Range informed by sorting success - lower than 97 %
HighAccuracy = find(Fraction_Match > 97);
Range_adjusted = [-5 -1:0.2:1.6 5]; 

% Preallocate - calculate fraction of components that are sorted acording
% to manual label assigned by Julian
Fraction_Match = nan(size(Range_adjusted, 2) - 1, 1);
test_output_score_sorted = output(idx, 1);
CompleteMatrix_sorted = CompleteMatrix(idx, 17);
test_score_sorted = test_output_score(idx, 1);

for x = 1:size(Range_adjusted, 2) - 1   
    IndexofRange = find(test_score_sorted(:, 1) > Range_adjusted(x) & ...
        test_score_sorted(:, 1) < Range_adjusted(x + 1));   
    Match = CompleteMatrix_sorted(IndexofRange) == 1 & test_output_score_sorted(IndexofRange) == 1 | ...
        CompleteMatrix_sorted(IndexofRange) == 0 & test_output_score_sorted(IndexofRange) == 0;
    Fraction_Match(x) = sum(Match)/numel(IndexofRange)*100;
end

figure; 
bar(Fraction_Match, 'k')
ylabel('Percent sorted corectly according to Julians labels')
xlabel('Distance from Decision Boundry')
box off
StringConcat = [cellfun(@num2str, mat2cell(Range_adjusted(1:end-1), 1, ones(numel...
    (Range_adjusted)-1, 1))', 'UniformOutput', false) repelem({'-'}, ...
    numel(Range_adjusted)-1)' cellfun(@num2str, mat2cell(Range_adjusted(2:end), 1, ...
    ones(numel(Range_adjusted)-1, 1))', 'UniformOutput', false)];
StringConcatFull = arrayfun(@(row)[(StringConcat{row,:})],(1:size(...
    StringConcat,1))','UniformOutput', false);
set(gca, 'XTick', 1:size(Range_adjusted, 2) - 1, 'XTickLabel', StringConcatFull)
xtickangle(45)

% Now sample the data in specified range
NumNeurons_perRange = round(Num_Sorting/(size(Range_adjusted, 2) - 1));
Inds_RandSample = cell(size(Range_adjusted, 2) - 1, 1);

for x = 1:size(Range_adjusted, 2) - 1   
    IndexofRange = find(test_output_score(:, 1) > Range_adjusted(x) & ...
        test_output_score(:, 1) < Range_adjusted(x + 1));  
    
    Inds_RandSample{x} = datasample(IndexofRange, NumNeurons_perRange, ...
        'Replace', false);
end

Inds_RandSample_concat = cat(1, Inds_RandSample{:});
Inds_RandSample_concat = Inds_RandSample_concat(randperm(length(...
    Inds_RandSample_concat)));

Dataset = CompleteNeurons(Inds_RandSample_concat(1:Num_Sorting), :);
Prediction = [output(Inds_RandSample_concat(1:Num_Sorting)) test_output_score(...
    Inds_RandSample_concat(1:Num_Sorting), 1)];
Indices_Selected = Inds_RandSample_concat(1:Num_Sorting);

% Next, save the Dataset
savefast(['\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Classification\TwoPhoton\Scoring_Results\Dataset_2P.mat'], 'Dataset');
savefast(['\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Classification\TwoPhoton\Scoring_Results\Prediction_2P.mat'], 'Prediction');
savefast(['\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Classification\TwoPhoton\Scoring_Results\Indices_2P.mat'], 'Indices_Selected');

clear CompleteNeurons
