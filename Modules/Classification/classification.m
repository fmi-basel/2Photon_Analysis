clear;
clc;

%% Load data
mat_data_train = load('Predictors.mat', 'Predictors_Tianlin').('Predictors_Tianlin');
mat_data_test = load('Predictors_Test.mat', 'Predictors_conc').('Predictors_conc');

X_train_raw = mat_data_train(:, 1:15);
X_test_raw = mat_data_test(:, 1:15);

X_train_raw(isinf(X_train_raw)) = 0;
X_test_raw(isinf(X_test_raw)) = 0;

% Label 0 indicates real neuron; 
y_train_raw = mat_data_train(:, end);
y_test_raw = mat_data_test(:, end);

test_num = length(y_test_raw);

%% Goal: Training a SVM classifier to maximize the precision. 
% The precision is defined as  # True Positive /(# True Positive + # False Positive), where
% # True Positive = number of cases where the classifier correctly predicted that a calcium trace IS NOT produced by a real neuron.
% # False Positive = number of cases where the classifier predicted that a calcium trace IS NOT from a real neuron, but in fact it IS from a real neuron. 

SVMModel = fitcsvm(X_train_raw,y_train_raw, 'Standardize',true,'KernelFunction','RBF','KernelScale','auto');
[train_output,train_output_score] = predict(SVMModel, X_train_raw);
precision_train = Evaluate_precision(y_train_raw, train_output);

[test_output,test_output_score] = predict(SVMModel, X_test_raw);
precision_test = Evaluate_precision(y_test_raw, test_output);
abs_test_score = abs(test_output_score(:, 1));


%% Discard some uncertain test data by checking their test_output_score; then re-evaluate the precision on the remaining part of the test data.

num_thres_points = 10; % the threshold for selecting uncertain data

thres_array = zeros(1,num_thres_points);
remaining_test_frac_array = zeros(1,num_thres_points);
test_precision_array = zeros(1,num_thres_points);

i = 0;

for thres = linspace(min(abs_test_score), max(abs_test_score) - 0.5, num_thres_points)
    i = i + 1;
    
    thres_indices = find(abs_test_score >= thres); % if the score is below the threshold, the data will be considered as uncertain and will be discarded.

    remaining_test_data_frac = length(thres_indices) / test_num; % the remaining fraction of test data after discarding the uncertain ones.

    test_precision = Evaluate_precision(y_test_raw(thres_indices), test_output(thres_indices));

    thres_array(i) = thres;
    remaining_test_frac_array(i) = remaining_test_data_frac;
    test_precision_array(i) = test_precision;
    
end

%%


subplot(2,1,1);
plot(thres_array, test_precision_array, 'LineWidth', 3)
grid on
xlabel('Data selection threshold',  'FontSize',14);
ylabel('Precision',  'FontSize',14);
ax = gca;
ax.FontSize = 14; 

subplot(2,1,2); 
plot(thres_array, remaining_test_frac_array, 'LineWidth', 3, 'Color', 'Black')
grid on
xlabel('Data selection threshold', 'FontSize',14);
ylabel({'Remaining fraction'; 'of test data'}, 'FontSize',14);

ax = gca;
ax.FontSize = 14; 

