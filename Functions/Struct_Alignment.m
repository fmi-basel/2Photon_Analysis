%% Align Structural Green and Red Channel
% For the Alignment the Reference Image is not good enough (at least not
% with GCaMP6f, need to use max Intensity projection

function [Out] = Struct_Alignment(FN_red, Max_proj)
    % How to do it:
    % First Record green + red simultanously
    % Now register green channel with the template of the motion correction
    % and apply the shifts to the red channel
    % Next step: Measure red intensity in the ROI, but I probbly need to
    % come up with a semiautomated way, since there might be cases where a
    % neuron has some red, but it is actually just a ROI that happens to be
    % close to a non labbelled GcAmp cell
    % Open 3 categories: excitatory, inhibitory and unsure
    % Bring up each cell and plot their fluorescence and ROI location + red
    
    
    
    Input_Red = mat2gray(loadtiff(FN_red));
    Movie_gray = mat2gray(Max_proj);
    Max_proj_calc = std(Movie_gray, [], 3);
    Max_proj_calc2 = max(Movie_gray, [], 3);
    figure;subplot(1,2,1);imshow(imadjust(Input_Red));subplot(1,2,2);imshow(imadjust(Max_proj_calc))    
    figure;subplot(1,2,1);imshow(imadjust(Max_proj_calc));subplot(1,2,2);imshow(imadjust(Max_proj_calc2))
    
    Cn = correlation_image(Movie_gray, 8, 512, 512, 0);
    
    saveastiff(Cn,'\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Results\2P\Julian\struct_Alignment\Test_green.tif');
    
    figure; subplot(1,3,1); imshow(imadjust(Max_proj_calc)); title('Standard deviation projection');
    subplot(1,3,2); imshow(Max_proj_calc2); title('Max Intensity projection'); 
    subplot(1,3,3); imshow(Cn); title('Local Pixel Correlation'); 
    
    T = adaptthresh(imadjust(Grayscale_red), 0.8, 'NeighborhoodSize', 5, 'Statistic', 'mean');CCbinary = find(T > 0.5);
    C = zeros(size(T, 1), size(T, 2)); C(CCbinary) = 1; L = bwlabel(C);

    figure;subplot(1,2,1);imagesc(T);colorbar;subplot(1,2,2);imagesc(L)
    Out = 0;

end


% 
% % % Some Random stuff for visualization
% Datta = reshape(Movie_gray, [512*512, size(Movie_gray, 3)]);
% Datta = reshape(movmean(Datta, 7, 2), [512 512 size(Movie_gray, 3)]);
% 
% figure;
% for kkk = 1:7:size(Movie_gray, 3)
%     imshow(Datta(:,:,kkk), [0 max(Datta(:))]);drawnow
% end


% 
% Moving_smoothed = mat2gray(imgaussfilt(A - imgaussfilt(A, 5), 2));
% figure;
% for i = 1:5:size(MC_Movies{1,1}, 3)
%     imagesc(MC_Movies{1,1}(:,:,i));drawnow
% end
%     
%          g = figure;
%         F(size(1:10:size(Moving_smoothed,3), 1)) = struct('cdata',[],'colormap',[]);
%         count = 0;
%         for ii = 1:10:size(Moving_smoothed,3)
%             count = count+1;
%             imshow(Moving_smoothed(:,:,ii), [0 max(Moving_smoothed(:))*0.4]);
%             drawnow;
%             F(count) = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%         end % caxis([0 max(Baseline_corrected(:))/4]);
%         
%         close(g);
%         v = VideoWriter('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Results\1P_Img\Example1');
%         open(v)
%         writeVideo(v, F)
%         close(v)   
    

