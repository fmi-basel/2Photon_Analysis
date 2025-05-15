%% Manual Validation - Scoring components

clear
close
clc

% Check if name is added
dims = [1 50];
Scorer_Name = inputdlg({''}, 'Please enter your Name ?', dims, {''});
Scorer_Name = lower(Scorer_Name{:});
Filepath = '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\shared\Miniscope_TwoPhotonscoring\';

% Load the dataset
load([Filepath 'Dataset_2P.mat'])

% First, check if there is already already some Selection
if ~exist([Filepath 'ScoringResults_2P\' Scorer_Name '_Scoring.mat'])
    ind_del = zeros(size(Dataset, 1), 1);
    m = 1;
    Scoring = struct;
else
    load([Filepath 'ScoringResults_2P\' Scorer_Name '_Scoring.mat'])
    m = Scoring.index_current;
    ind_del = Scoring.Deletion;
end

User_termination = false;
while m <= size(Dataset, 1) && ~User_termination
    
    % time
    T = size(Dataset{m, 2}, 2);
    t = 1:T;
    t = t/(Dataset{m, 6});
    str_xlabel = 'Time (Sec.)';
    
    % Plot the component that enables nice selection
    subplot(2, 4, [1 5]); cla;
    imagesc(Dataset{m, 1}); colormap gray;hold on;
    set(gca, 'XTick', [], 'YTick', []) 
    
    % temporal components
    subplot(2,4, 2:4); cla;
    plot(t, movmean(zscore(Dataset{m, 3}), ceil(Dataset{m, 6}*2)), 'k', 'linewidth', 0.5); hold on;
    ylim([min(zscore(Dataset{m, 3})), ...
        max(zscore(Dataset{m, 3}))]);
    box off; ylabel('Fluorescence'); set(gca, 'XTick', [])
    
    try
        vline(600:600:size(Dataset{m, 3}, 2)/(Dataset{m, 6}), 'b', [], 1);  
    catch
        vline(600:600:600, 'b', [], 1); 
    end
    
    xlim([t(1), t(end)]);
    if ind_del(m)
        title('Currently deleted', 'color', 'r'); 
    else
        title('Currently not deleted');
    end

    suptitle(sprintf('Neuron %d out of %d', m, size(Dataset, 1)));
    
    % temporal components
    subplot(2, 4, 6:8);cla;
    plot(t, zscore(Dataset{m, 3}), 'k', 'linewidth', 0.5); hold on;
    plot(t, Dataset{m, 2}, 'r', 'linewidth', 1); hold on;
     ylim([min(zscore(Dataset{m, 3})), ...
        max(zscore(Dataset{m, 3}))]);
    xlabel(str_xlabel); box off; ylabel('Fluorescence')
    try
        vline(600:600:size(Dataset{m, 3}, 2)/(Dataset{m, 6}), 'b', [], 1);  
    catch
        vline(600:600:600, 'b', [], 1); 
    end
    xlim([t(1), t(end)]);


    % Now do selection
    fprintf('Neuron %d, keep(k, default)/delete(d)/backward(b)/Terminate Component Selection(e):    ', m);
    temp = input('', 's');
    
    if temp=='d'
        ind_del(m) = true;
        m = m+1;
    elseif strcmpi(temp, 'b')
        if m == 1
            m = 1;
        else
            m = m - 1;
        end        
    elseif strcmpi(temp, 'k') 
        ind_del(m) = false;
        m = m+1;
    elseif strcmpi(temp, 'e')
        User_termination = true;
    else
        m = m + 1;
    end

end

Scoring.index_current = m;
Scoring.Deletion = ind_del;

% Now save the Output
save([Filepath 'ScoringResults_2P\' Scorer_Name '_Scoring.mat'], 'Scoring');
