%% ROI Segmentation

function [obj] = Corr_Img_Segmentation(A, vis, num_planes)

    % Calculate Correlation Image
    % A is data in format X,Y,T
    
    DS_Dat_vis = [];
    load([A '\processed_data_Folder\DS_Dat_vis.mat']);
    Sizes = size(DS_Dat_vis);
    
    for i = 1:num_planes
        if num_planes == 1
            Cn = correlation_image(mat2gray(DS_Dat_vis), 8, Sizes(1), Sizes(2), 0);
            % Max_proj = max(A(:,:,1:10:end), [], 3);
            STD_Proj = std(mat2gray(DS_Dat_vis), [], 3);
            Median_Img = median(DS_Dat_vis, 3);
        else
            Cn = correlation_image(mat2gray(DS_Dat_vis(:,:,:, i)), 8, Sizes(1), Sizes(2), 0);
            % Max_proj = max(A(:,:,1:10:end), [], 3);
            STD_Proj = std(mat2gray(DS_Dat_vis(:,:,:,i)), [], 3);
            Median_Img = median(DS_Dat_vis(:,:,:,i), 3);
        end
        Cn(isnan(Cn)) = 0;
        Merge = STD_Proj + mat2gray(Median_Img) + 2 * Cn;

        if vis
            figure;
            subplot(1,4,1); imagesc(Cn); axis off; subplot(1,4,2); imagesc(STD_Proj); axis off
            subplot(1,4,3); imagesc(Median_Img); axis off;subplot(1,4,4); imagesc(Merge); axis off
        end

        % save relevant parameters
        obj.(['Cor_Images_' num2str(i)]).Corr_Img = cat(3, Cn, STD_Proj, mat2gray(Median_Img), Merge);
    end
end