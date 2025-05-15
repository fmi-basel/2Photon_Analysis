%% Frame Timings

function [Timing, FS] = FrameTS(Files, num_planes, Scanning_mode, Length_Recording)
    
    % Save Timings of Images
    searchString = fullfile(Files, '*.xml');
    Filenames = dir(searchString);
    
    index_Triggers = ~cellfun(@(x) any(strfind(lower(x), 'voltage')),{Filenames(:).name}) & ...
           transpose(cat(1, Filenames(:).bytes) > 2000);
    
    if ~any(index_Triggers)
        index_Triggers = ~cellfun(@(x) any(strfind(lower(x), 'voltage')),{Filenames(:).name});
    end
    
    [s] = xml2struct2([Files Filenames(index_Triggers).name]);      
    format long
    
    if num_planes > 1
        Sequences_Img = cellfun(@(c) {c.Frame}, s.PVScan.Sequence, 'UniformOutput', true);
        
        Timing_abs = nan(Length_Recording, num_planes);
        Timing_rel = nan(Length_Recording, num_planes);
        
        for i = 1:num_planes
            Timing_rel(:,i) = cellfun(@str2double, cellfun(@(c) {c{i}.Attributes.relativeTime}, Sequences_Img(1:Length_Recording), 'UniformOutput', true));
            Timing_abs(:,i) = cellfun(@str2double, cellfun(@(c) {c{i}.Attributes.absoluteTime}, Sequences_Img(1:Length_Recording), 'UniformOutput', true));
        end
        

        switch Scanning_mode
            case 'Bidirectional'
                Timing_rel(2:2:end, :) = fliplr(Timing_rel(2:2:end, :));
                Timing_abs(2:2:end, :) = fliplr(Timing_abs(2:2:end, :));
%                 Timing_abs(2:2:end, :) = fliplr(Timing_rel(2:2:end, :));
                
                Timing = permute(cat(3, Timing_rel, Timing_abs), [1, 3, 2]);
            case 'SingleDirection'
                Timing = permute(cat(3, Timing_rel, Timing_abs), [1, 3, 2]); 
        end
       
    else
        Sequences_Img = s.PVScan.Sequence.Frame;
        
        Timing = transpose([cell2mat(cellfun(@(c) str2double(c.Attributes.relativeTime), Sequences_Img, 'UniformOutput', false)); ... 
             cell2mat(cellfun(@(c) str2double(c.Attributes.absoluteTime), Sequences_Img, 'UniformOutput', false))]);

    end
    
    FS = round(size(Sequences_Img, 2)/Timing(end, 1, num_planes));

end

% Old way that failed, if the Voltage Output is bigger than the frame trigger 
%     [~, idx_time] = max(cat(1, Filenames(:).bytes));
%     [s] = xml2struct2([Files Filenames(idx_time).name]);

