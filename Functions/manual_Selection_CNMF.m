function [ind_del_sel, ind_del, Data] = manual_Selection_CNMF(Data, FS, Corr_Imgs, RegionSub)


%% view all components and delete components manually. it shows spatial
%   components in the full-frame and zoomed-in view. It also shows temporal
%   components
%% input:
%   ind: vector, indices of components to be displayed, no bigger than the maximum
%       number of neurons
%   C2:  K*T matrix, another temporal component to be displayed together
%       with the esitmated C. usually it is C without deconvolution.

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016
% Completely remodelled Julian Hinz, FMI 2018


ind_del = false(size(Data.C, 1), 1);     % indicator of deleting neurons
gSiz = 8;        % maximum size of a neuron
Cn = Data.Cn;
ctr = Data.Centroids;      %neuron's center

% time
T = size(Data.C, 2);
t = 1:T;
if ~isnan(FS)
    t = t/FS;
    str_xlabel = 'Time (Sec.)';
else
    str_xlabel = 'Frame';
end

To_exclude_weird = logical(zeros(1, size(Data.Boundries, 1)));
g = figure;
imagesc(Corr_Imgs); 
colormap gray;hold on;
for mm = 1:size(Data.Boundries, 1)
    try
        cont = Data.Boundries{mm}; 
        plot(cont{1}(:,1),cont{1}(:,2),'Color', 'b', 'linewidth', 0.75);
    catch
        To_exclude_weird(mm) = true;
    end
end

F = getframe ;
[X,~] = frame2im(F);
All_comp = imresize(X,[size(Corr_Imgs)], 'bicubic');
close(g)

if RegionSub
    Mean_Calc = mean(Corr_Imgs(:));
    g = figure('position', [1000, 300, 800, 800]);
    title('Please select Region/Cells to be included by drawing a ROI around them !'); 
    imagesc(Corr_Imgs, [0 Mean_Calc + 0.3*Mean_Calc]); 
    colormap gray;hold on;

    for mm = 1:size(Data.Boundries, 1)
        if ~To_exclude_weird(mm)
            cont = Data.Boundries{mm}; 
            plot(cont{1}(:,1),cont{1}(:,2),'Color', 'b', 'linewidth', 0.75);
        end
    end
    set(gca, 'XTick', [], 'YTick', [])
    title('Please select Region/Cells to be included by drawing a ROI around them !'); 
    Region_selected = roipoly();
    Region_selected_cell = double(repmat(reshape(Region_selected, size(Region_selected, 1) * ...
        size(Region_selected, 2), 1), 1, size(Data.A, 2)));
    Region_circled = (sum(Region_selected_cell - double(Data.A > 0) < 0, 1) ./ sum(Data.A > 0, 1)) > 0.2;

    To_exclude = Region_circled | To_exclude_weird;

    close(g)
else
    To_exclude = To_exclude_weird;
end

%% start viewing neurons
g = figure('position', [50, 390, 2000, 600]);
ind_del(To_exclude) = true;
Indices_delete = find(~To_exclude);

m=1;
while and(m>=1, m<=size(Indices_delete, 2))
    
    %% Position in Image
    subplot(341); cla;
    cont = Data.Boundries{Indices_delete(m)}; 
    imagesc(Corr_Imgs); colormap gray;hold on;
    plot(cont{1}(:,1),cont{1}(:,2),'Color', 'r', 'linewidth', 1.5); 
    axis equal; axis off;
    hold off
    
    
    %% Position in Image
    subplot(342); cla;
    cont = Data.Boundries{Indices_delete(m)}; 
    imagesc(Corr_Imgs); colormap gray;hold on;
    plot(cont{1}(:,1),cont{1}(:,2),'Color', 'r', 'linewidth', 1.5); 
    axis equal; axis off;
    hold off
    
    x0 = ctr(Indices_delete(m), 2);
    y0 = ctr(Indices_delete(m), 1);
    
    if ~isnan(x0)
        xlim(x0+[-gSiz, gSiz]*2);
        ylim(y0+[-gSiz, gSiz]*2);
    end

    %% Position in Image next to all others
    subplot(343); cla;
    imagesc(All_comp); hold on;
       
    plot(cont{1}(:,1),cont{1}(:,2),'Color', 'r', 'linewidth', 1.5); 
    hold off
    axis equal; axis off;
    
    %% Position in Image next to all others
    subplot(344); cla;
    imagesc(All_comp); hold on;
       
    plot(cont{1}(:,1),cont{1}(:,2),'Color', 'r', 'linewidth', 1.5); 
    hold off
    axis equal; axis off;
    
    if ~isnan(x0)
        xlim(x0+[-gSiz, gSiz]*2);
        ylim(y0+[-gSiz, gSiz]*2);
    end

    
    %% temporal components
    subplot(3,4,5:8);cla;
    plot(t, movmean(zscore(Data.C(Indices_delete(m), :)), ceil(FS*2)), 'k', 'linewidth', 0.5); hold on;
    ylim([min(zscore(Data.C(Indices_delete(m), :))), ...
        max(zscore(Data.C(Indices_delete(m), :)))]);
     box off; ylabel('Fluorescence'); set(gca, 'XTick', [])
    
    if ind_del(Indices_delete(m))
        title('Currently deleted', 'color', 'r'); 
    else
        title('Currently not deleted');
    end
    vline(600:600:size(Data.C, 2)/FS, 'r', [], 3); 
    xlim([t(1), t(end)]);
    suptitle(sprintf('Neuron %d out of %d', Indices_delete(m), size(Data.C, 1)));
    
    %% temporal components
    subplot(3,4,9:12);cla;
    plot(t, zscore(Data.C(Indices_delete(m), :)), 'k', 'linewidth', 0.5); hold on;
     ylim([min(zscore(Data.C(Indices_delete(m), :))), ...
        max(zscore(Data.C(Indices_delete(m), :)))]);
    xlabel(str_xlabel); box off; ylabel('Fluorescence')
    vline(600:600:size(Data.C, 2)/FS, 'r', [], 3); 
    xlim([t(1), t(end)]);
    %% save images
    fprintf('Neuron %d, keep(k, default)/delete(d)/cancel(tc)/delete all(da)/backward(b)/end(e):    ', m);
    temp = input('', 's');
    
    if temp=='d'
        ind_del(Indices_delete(m)) = true;
        m = m+1;
    elseif strcmpi(temp, 'b')
        if m == 1
            m = 1;
        else
            m = m - 1;
        end        
    elseif strcmpi(temp, 'da')
        ind_del(Indices_delete(m):end) = true;
        break;
    elseif strcmpi(temp, 'k')
        ind_del(Indices_delete(m)) = false;
        m= m+1;
    elseif strcmpi(temp, 'e')
        break;
    else
        m = m + 1;
    end
    
end
close(g)

ind_del_sel = ind_del;

%%% Now delete Cells that are close and highly correlated
% First, check if there were already cells deleted and delete the
% corresponding matches
Overlapping_del_Corr = intersect(unique(Data.Merge_DelInd_Sub), find(ind_del));
ToDisplay = Data.Merge_DelInd_Sub(~any(ismember(Data.Merge_DelInd_Sub, ...
    Overlapping_del_Corr), 2), :);

% Check the 3 scenarios, there are no cells left with overlap, all the
% cells are unique, there are cells that are close and correlated to
% several other cells

if isempty(ToDisplay)
    return
else
    % First, identify how many cells that are overlapping with more than
    % one other cell - solve the problem if there is a cell that is
    % overlapping with 2 other cells
    
    [Identity, ~, frequent_cellnum] = unique(ToDisplay, 'stable');
    Frequency = accumarray(frequent_cellnum, 1);

    New_Corr_dist = nan(size(ToDisplay, 1), max(Frequency) + 1);
    z = 1;

    while ~isempty(Identity)
       Corr_temps = unique(ToDisplay(any(ToDisplay == Identity(1), 2), :));
       New_Corr_dist(z, 1:numel(Corr_temps)) = Corr_temps;
       Identity(ismember(Identity, Corr_temps)) = []; 
       z = z + 1;
    end

    New_Corr_dist(all(isnan(New_Corr_dist), 2), :) = [];
    Colours = distinguishable_colors(size(New_Corr_dist, 2));
    
    warning('off','all');
    
    g = figure('position', [50, 390, 1000, 600]);
    for z = 1:size(New_Corr_dist, 1)
        Temp = zeros(size(Cn, 1)*size(Cn, 2), 1);
        for zz = 1:sum(~isnan(New_Corr_dist(z, :)))
            Temp(Data.A(:, New_Corr_dist(z, zz)) > 0) = zz;            
        end
        cc = label2rgb(reshape(Temp, size(Cn, 1), size(Cn, 2)), ...
            Colours(1:sum(~isnan(New_Corr_dist(z, :))), :));
        subplot(121);
        imagesc(Cn); colormap bone; hold on;imagesc(cc, 'AlphaData', .5); box off;
        set(gca, 'XTick', [], 'YTick', []);
        subplot(122)
        h = cell(sum(~isnan(New_Corr_dist(z, :))), 1);
        hold on;
        for zz = 1:sum(~isnan(New_Corr_dist(z, :)))
            bb.(['Field' num2str(zz)]) = plot(t, Data.C( New_Corr_dist(z, zz), :), 'color', Colours(zz, :));
        end
        xlabel(str_xlabel); ylabel('Fluorescence'); box off
        legend('Cell 1', 'Cell 2', 'Cell 3'); xlim([0 max(t) + 5])
        drawnow
        
        fprintf('Keep all neurons (k, default), Delete all neurons (d):    ');
        temp = input('', 's');
        
        for zz = 1:sum(~isnan(New_Corr_dist(z, :)))
            delete(bb.(['Field' num2str(zz)]));
        end
        
        if isempty(temp) || strcmpi(temp, 'k')
            continue   
        elseif strcmpi(temp, 'd')
            ind_del(New_Corr_dist(z, ~isnan(New_Corr_dist(z, :)))) = 1;
        else
            try
                Numbers = cellfun(@str2num, regexp(temp, ['[1-' ...
                              num2str(sum(~isnan(New_Corr_dist(z, :)))) ']'], 'match'));
                ind_del(New_Corr_dist(z, Numbers)) = 1;
            catch
                disp('Selection failed, both neurons are still preserved !')
            end
        end
        
    end
    
    warning('on','all');
    
    close(g)
    return;

    
end




end

% Old Plot
% 
% 
%    %% full-frame view
%     subplot(241); cla;
%     % obj.image(obj.A(:, ind(m))); %
%     imagesc(reshape(Data.A(:, Indices_delete(m)), [size(Cn, 1), size(Cn, 2)]));
%     axis equal; axis off; colormap gray
% 
%     %% zoomed-in view
%     subplot(242); cla;
%     imagesc(reshape(Data.A(:, Indices_delete(m)), [size(Cn, 1), size(Cn, 2)])); %
%     axis off;axis equal;
% 
%     x0 = ctr(Indices_delete(m), 2);
%     y0 = ctr(Indices_delete(m), 1);
%     
%     if ~isnan(x0)
%         xlim(x0+[-gSiz, gSiz]*2);
%         ylim(y0+[-gSiz, gSiz]*2);
%     end
%     
%     if ind_del(Indices_delete(m))
%         title(sprintf('Neuron %d out of %d', Indices_delete(m), size(Data.C, 1)), 'color', 'r');
%     else
%         title(sprintf('Neuron %d out of %d', Indices_delete(m), size(Data.C, 1)));
%     end
%     
%     %% Position in Image
%     subplot(243); cla;
%     cont = Data.Boundries{Indices_delete(m)}; 
%     imagesc(Corr_Imgs); colormap gray;hold on;
%     plot(cont{1}(:,1),cont{1}(:,2),'Color', 'r', 'linewidth', 1.5); 
%     axis equal; axis off;
%     hold off
%     
%     %% Position in Image next to all others
%     subplot(244); cla;
%     imagesc(All_comp); hold on;
%        
%     plot(cont{1}(:,1),cont{1}(:,2),'Color', 'r', 'linewidth', 1.5); 
%     hold off
%     axis equal; axis off;
%     
%     %% temporal components
%     subplot(2,4,5:8);cla;
%     plot(t, zscore(Data.C(Indices_delete(m), :)), 'k', 'linewidth', 0.5); hold on;
%     xlim([t(1), t(end)]); ylim([min(zscore(Data.C(Indices_delete(m), :))), ...
%         max(zscore(Data.C(Indices_delete(m), :)))]);
%     xlabel(str_xlabel); box off; ylabel('Fluorescence')




