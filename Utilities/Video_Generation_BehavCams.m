%% Generate Video that contains LickSequences from a few example animals

% Clean slate
clear all
close all
clc

addpath(genpath('\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Github\gluthi-2Photon'))

% Video Location
Vid_Location = '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\shared\JulianDatatransfer\Videos\';
Vid_Name = 'Camera_2';

% Excel Location
Behav_Location = '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\shared\JulianDatatransfer\Bruker_Recording\';

% Sampling frequency Bruker
Sampling_Freq_V = 5000;
VideoFR = 194.1;

% Now to processing 
VideoFiles = dir(fullfile(Vid_Location, '*.mp4'));
VideoNames = {VideoFiles.name};

AnimalIdent = cell2mat(cellfun(@(c) str2double(c(3:9)), VideoNames, 'UniformOutput', false))';
Animals_un = unique(AnimalIdent);

Video_frames_Animal = cell(numel(Animals_un), 1);

for x = 1:numel(Animals_un)
    
    SubselectedVids = VideoNames(AnimalIdent == Animals_un(x));
    Sub_Rightcam = SubselectedVids(cell2mat(cellfun(@(c) contains(c, Vid_Name), ...
        SubselectedVids, 'UniformOutput', false)));
    Video_Frames_Ses = cell(numel(Sub_Rightcam), 4);
    
    for xx = 1:numel(Sub_Rightcam)
        tmp_vid = Sub_Rightcam{xx};
        
        % First identify which session and extract Licking related info
        find_underscore = strfind(tmp_vid, '_');
        Ses_ID = str2double(tmp_vid(find_underscore(3)+1:find_underscore(4)-1));
        
        if Ses_ID <= 2
            continue
        end
        
        Location = [Behav_Location 'A_' num2str(Animals_un(x)) '\Session_' ...
            num2str(Ses_ID) '\Imaging\']; 
        
        % Now proceed to find Folder
        Inside_Fold = dir(Location);
        Inside_Fold = {Inside_Fold.name};
        
        try
            Inside_Fold_ext = Inside_Fold{cell2mat(cellfun(@(c) size(c,2),  Inside_Fold, ...
                'UniformOutput', false)) > 10};
        catch
           disp('couldnt identify only one Folder, please check if there is multiple Folders') 
           continue
        end
        
        % Now find the excelsheet and process it
        Location_full = [Location Inside_Fold_ext '\'];
        try
            [Timing_maxLick] = min_Analze_Trigger_File(Location_full, Sampling_Freq_V);
        catch
           continue 
        end

        
        if numel(Timing_maxLick) < 3
           continue 
        end
        
        Video_Frames_Ses(xx, 2:end) = transpose(cellfun(@(c) [round((mean(c)-10)*VideoFR) round((mean(c)...
            +10)*VideoFR)], Timing_maxLick, 'UniformOutput', false));
        Video_Frames_Ses{xx, 1} = [Vid_Location tmp_vid];
        disp(['Finished ' num2str(xx)])
    end
    
    Video_Frames_Ses(all(cellfun(@isempty, Video_Frames_Ses), 2), :) = [];
    Video_frames_Animal{x} = Video_Frames_Ses;
    disp(['Finished Animal ' num2str(x)])

end

Combined_Mat = cat(1, Video_frames_Animal{:});
% Convert cell to a table and use first row as variable names
T = cell2table(Combined_Mat);
 
% Write the table to a CSV file
writetable(T,'\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\shared\JulianDatatransfer\LickBouts.csv')








