%% Read and analyze the Trigger file
% Might only work in cases with short trace, need to check

function [Stimulation, Licks, RunningSpeed, Pump_activity, Camera_TTL] = ...
    Analze_Trigger_File(Input_Fold, Sampling_Freq_V, tick2cm, Video_Sampling, Duty_cycle)
    
    % Find the .csv file in specified Folder
    Temporary = dir([Input_Fold '\Imaging\']);
    Temporary = Temporary(cellfun(@numel, {Temporary.name}) > 10);
    
    if size(Temporary, 1) == 1
        Files_Img = [Input_Fold '\Imaging\' Temporary.name '\'];
    elseif size(Temporary, 1) > 1 && size(Temporary, 1) < 1000 
        error('Check whats weird here !')
    else
        Files_Img = [Input_Fold '\Imaging\']; 
    end
    
    Voltage_Recording = fullfile(Files_Img, '*.csv');
    Filenames = dir(Voltage_Recording);
   
    M = csvread([Filenames.folder '\' Filenames.name], 1, 0);
    opts = detectImportOptions([Filenames.folder '\' Filenames.name]);
    Header = opts.VariableNames; 
    
    % xlsread caused problems, therefore changed
    % [~,Header,~] = xlsread([Filenames.folder '\' Filenames.name], 1) ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Triggers for events
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Triggerf = cellfun(@(x) any(strfind(lower(x), 'trig')), Header);
    
    if any(Triggerf)
        Indices = find(diff(M(:, Triggerf)) > 0.1 | diff(M(:, Triggerf)) < -0.1);
        Events = ~round(M(Indices + 1, Triggerf), 1, 'decimals') == 0;
        
        cl_Ind = Indices(Events);
        Start_Stim = cl_Ind/Sampling_Freq_V;
        
        Length_Stim = diff(Indices)/Sampling_Freq_V;
        Length_Stim = Length_Stim(Events(1:end-1));

        % Length_Stim = transpose(diff(reshape(Indices, [2, size(Indices, 1)/2])))/Sampling_Freq_V;
        Stim_Identity = round(M(cl_Ind + 1, Triggerf), 1, 'decimals');
        
        if numel(Start_Stim) ~= numel(Length_Stim)
            Stimulation = [Start_Stim(1:numel(Length_Stim)) Length_Stim Stim_Identity(1:numel(Length_Stim))];
        else
            Stimulation = [Start_Stim Length_Stim Stim_Identity];
        end
    else
        Stimulation = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Licks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % check if Licks were recorded
    Licksind = cellfun(@(x) any(strfind(x, 'Licks')), Header);
    
    if any(Licksind)
        
        % find on and offset of Licks
        Licks_calc_on = find(diff(M(:, Licksind) < - 0.5) == 1) + 1;
        Licks_calc_off = find(diff(M(:, Licksind) > 0.5) == 1) + 1;
        
        if size(Licks_calc_on, 1) > 10


            % Now we need to discard of artifacts
            Licks_on_off = nan(size(Licks_calc_on, 1), 2);
            for z = 1:size(Licks_calc_on, 1)

                clos_match = Licks_calc_off - Licks_calc_on(z) < 500 ...
                    & Licks_calc_off - Licks_calc_on(z) > 0;

                if sum(clos_match) == 1
                    Licks_on_off(z, :) = [Licks_calc_on(z) Licks_calc_off(clos_match)];
                elseif sum(clos_match) > 1
                    idx = find(clos_match);
                    if Licks_calc_off(idx(1)) -  Licks_calc_on(z) < 100
                        Licks_on_off(z, :) = [Licks_calc_on(z) Licks_calc_off(idx(1))];
                    end
                end
            end

            Licks_on_off(isnan(Licks_on_off(:,1)), :) = [];

            Lick_ON_OFF = [Licks_on_off(:,1) Licks_on_off(:,2) ...
                Licks_on_off(:,2) - Licks_on_off(:,1)]/Sampling_Freq_V; % Calculate Onset, Offset and Duration    

    %         figure; hist(Licks_calc_on2 - Licks_calc_on);
    %         abs(M(Licks_calc_on, Licksind));
    %         abs(M(Licks_calc_off, Licksind));
    %         
    %         Licks_calc_on2 = find(diff(M(:, Licksind) < - 0.3) == -1) + 1;
    %         Licks_calc_off2 = find(diff(M(:, Licksind) > 0.3) == -1);
    %         figure; hist(Licks_calc_off2 - Licks_calc_off)
    %         
    %         Vals = M(Licks_calc_on, Licksind);
    %         Takeout = Vals < mean(Vals) - 2*std(Vals);
    %         sum(Takeout)
    %         
    %         Vals2 = M(Licks_calc_off, Licksind);     
    %         Takeout2 = Vals2 < mean(Vals2) - 2*std(Vals2);
    %         
    %         sum(Takeout2)
    %         figure; hist(abs(Licks_calc_on(~Takeout) - Licks_calc_off), [0:1:5000]);
    %         
%             figure;plot(M(:, Licksind), 'k');
%             hold on;
%             scatter(Licks_on_off(:,2), zeros(size(Licks_on_off, 1), 1) +0.3, 'r')
%             scatter(Licks_on_off(:,1), zeros(size(Licks_on_off, 1), 1) - 0.3, 'b')
%             title([Temporary.name], 'Interpreter', 'none')
    %         
    %         % Calculate Duration and error out if there is an uneven number of
    %         % on and offsets
    %         if numel(Licks_calc_on) == numel(Licks_calc_off)
    %             Lick_ON_OFF = [Licks_calc_on Licks_calc_off ...
    %                 Licks_calc_off - Licks_calc_on]/Sampling_Freq_V; % Calculate Onset, Offset and Duration            
    %         elseif numel(Licks_calc_on) > numel(Licks_calc_off)
    %             figure; hist(abs(Licks_calc_on(1:end-1) - Licks_calc_off), [0:1:20000]);
    %             test = find(abs(Licks_calc_on(1:end-1) - Licks_calc_off) > 150);
    %             
    %             Licks_calc_on(test(1)) = [];
    %             figure;hist(Licks_calc_on- Licks_calc_off, [-500:1:20000]);
    %         else
    %             error('Check this again !')
    %         end
    %         
            % Fit normal distribution across the Lick duration to exclude Licks
            % that seem to be flawed
            pd = fitdist(Lick_ON_OFF(:,3), 'Normal');
            Exclude_Licks_Dur = Lick_ON_OFF(:,3) <= pd.mu - 4*pd.sigma | ...
                Lick_ON_OFF(:,3) >= pd.mu + 4*pd.sigma; % exclude licks at the edge of the distribution
            idx_ILI = [0; diff(Lick_ON_OFF(:, 1)) < 0.08]; % exclude licks with an ITI of less than 0.8 s 

            exclusion = Exclude_Licks_Dur + logical(idx_ILI);
            exclusion(exclusion > 1) = 1; % combine the results

            Cl_Lick_ON_OFF = Lick_ON_OFF(~logical(exclusion), :); 

            % Plot the resulting distribution and the excluded Licks

            Norm_Gaus = max(hist(Lick_ON_OFF(:, 3), [0:0.001:0.2])) / ...
                max(normpdf(0:0.001:0.2, pd.mu, pd.sigma));

            g = figure('position', [500 100 1500 1000]); hist(Lick_ON_OFF(~logical(exclusion), 3),[0:0.001:0.2])
            hold on; hist(Lick_ON_OFF(logical(exclusion), 3), [0:0.001:0.2])
            plot(0:0.001:0.2, normpdf(0:0.001:0.2, pd.mu, pd.sigma)*Norm_Gaus, 'k', 'linewidth', 1.5)
            h = findobj(gca,'Type','patch');
            h(1).FaceColor = 'r'; h(1).EdgeColor = 'w';
            h(2).FaceColor = 'b'; h(2).EdgeColor = 'w';
            vline([pd.mu - 2*pd.sigma; pd.mu + 2*pd.sigma], 'k-', '', 1.5)
            xlabel('Lick Duration [s]'); ylabel('# Frequency'); box off

            saveas(g, [Input_Fold '\processed_data_Folder\Lick_Duration_Distribution.jpg']);
            close(g)        
            
            figure;plot([1:size(M, 1)]/Sampling_Freq_V, M(:, Licksind), 'k');
            hold on;
            scatter(Cl_Lick_ON_OFF(:,2), zeros(size(Cl_Lick_ON_OFF, 1), 1) + 0.5, 'r')
            scatter(Cl_Lick_ON_OFF(:,1), zeros(size(Cl_Lick_ON_OFF, 1), 1) - 0.5, 'b')
            title([Temporary.name], 'Interpreter', 'none')

            
%             figure
%             ax1 = subplot(1,2,1);
%             plot([1:size(M, 1)]/Sampling_Freq_V, M(:, Licksind), 'k');
%             hold on;
%             scatter(Licks_on_off(:,2)/Sampling_Freq_V, zeros(size(Licks_on_off, 1), 1) +0.3, 'r')
%             scatter(Licks_on_off(:,1)/Sampling_Freq_V, zeros(size(Licks_on_off, 1), 1) - 0.3, 'b')
%             title([Temporary.name], 'Interpreter', 'none')
%             xlim([0 3000])
%             ax2 = subplot(1,2,2);
%             plot([1:size(M, 1)]/Sampling_Freq_V, M(:, Licksind), 'k');
%             hold on;
%             scatter(Cl_Lick_ON_OFF(:,2), zeros(size(Cl_Lick_ON_OFF, 1), 1) + 0.5, 'r')
%             scatter(Cl_Lick_ON_OFF(:,1), zeros(size(Cl_Lick_ON_OFF, 1), 1) - 0.5, 'b')
%             title([Temporary.name], 'Interpreter', 'none')
%             xlim([0 3000])
%             linkaxes([ax1,ax2],'x');

            % Save the Lick related data 
            Licks.Cl_Lick_ON_OFF = Cl_Lick_ON_OFF;
            Licks.Lick_ON_OFF = Lick_ON_OFF;
            Licks.Licks_calc_on = Licks_calc_on/Sampling_Freq_V;
            Licks.Licks_calc_on = Licks_calc_off/Sampling_Freq_V;

            % Now go on to identify Lick bouts
            Frequency_Licks = diff(Cl_Lick_ON_OFF(:, 1));
            g = figure; hold on;hist(Frequency_Licks(Frequency_Licks < 0.3), [0:0.005:0.3]);
            xlim([0 0.35]); vline([median(Frequency_Licks(Frequency_Licks < 0.3))], 'r-', '', 1.5)
            h = findobj(gca,'Type','patch');
            h(1).FaceColor = 'k'; h(1).EdgeColor = 'w';
            xlabel('Inter Lick Interval'); ylabel('# Frequency')
            saveas(g, [Input_Fold '\processed_data_Folder\Lick_ITI_Distribution.jpg']);
            close(g)

            Licks.MedFreq = 1/median(Frequency_Licks); % in Hz

            % Bout Analysis       
            [Bout_on, Bout_offset] = find_sequences_Licks(diff([0; Cl_Lick_ON_OFF(:, 1)]) < 0.3, 1);

            if ~isempty(Bout_on)
                Bout_on = Bout_on - 1;

                if Bout_on(1) == 0
                    Bout_on(1) = 1;
                end

                Bout_Collection = cell(size(Bout_on, 2), 1);
                Median_ILI = nan(size(Bout_on, 2), 1);
                Mean_frequency = nan(size(Bout_on, 2), 1);
                Bout_length = (Bout_offset - Bout_on) + 1;

                for z = 1:size(Bout_on, 2)
                    Bout_Collection{z} = Cl_Lick_ON_OFF(Bout_on(z):Bout_offset(z), 1);
                    Median_ILI(z) = median(diff(Bout_Collection{z}));
                    Mean_frequency(z) = (Cl_Lick_ON_OFF(Bout_offset(z)) - Cl_Lick_ON_OFF(Bout_on(z)))/numel(Bout_on(z):Bout_offset(z))  ;
                end
            else
                Bout_Collection = nan;
                Median_ILI = nan;
                Bout_length = nan;
                Mean_frequency = nan;
                Bout_on = nan;
                Bout_offset = nan;
            end

            Licks.Bout_Collection = Bout_Collection;
            Licks.Median_ILI = Median_ILI;
            Licks.Bout_length = Bout_length;
            Licks.Mean_frequency = Mean_frequency;    
            Licks.Bout_On_Off = [Bout_on; Bout_offset]';
        else
            Licks = [];
        end
        
    else
        Licks = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Running Speeed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    RS_Rec = cellfun(@(x) any(strfind(lower(x), 'rec')), Header) | ...
        cellfun(@(x) any(strfind(lower(x), 'run')), Header) | ...
        cellfun(@(x) any(strfind(lower(x), 'speed')), Header);
       
    if any(RS_Rec)
        Binarized_Transitions = M(:, RS_Rec);
        Indices = round(Binarized_Transitions) > 0;
        Binarized_Transitions(Indices) = 1;
        Binarized_Transitions(~Indices) = 0;

        % Lastly calculate a moving sum over the Transitions and convert to cm
        Transitions = [0; diff(Binarized_Transitions,1,1) ~= 0];
        RunningSpeed = movsum(Transitions, Sampling_Freq_V) / tick2cm;
    else
        RunningSpeed = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Camera TTLs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Cam = cellfun(@(x) any(strfind(lower(x), 'cam')), Header);
    
    if any(Cam)
        if sum(Cam) > 1
            % Find Position in Bruker Recording and sort by number of
            % camera
            CamPos = find(Cam);
            Header_sub = Header(CamPos);
            PosShuff = nan(numel(CamPos), 1);
            for j = 1:numel(CamPos)
                PosShuff(j) = CamPos(cell2mat(cellfun(@(c) any(strfind(['camera_' ...
                    num2str(j)], lower(c))), Header_sub, 'UniformOutput', false)));
            end
            
            % Run through camera and identify the TTLs
            Camera_TTL = cell(numel(CamPos), 1);
            for j = 1:numel(CamPos)
                
                % Identify 0 crossing of differential pos/neg
                CamTTL = round(M(:, PosShuff(j)));
                TTL_High = (find(diff(CamTTL,1,1) > 0.5) + 1)/Sampling_Freq_V;
                TTL_Low = (find(diff(CamTTL,1,1) < -0.5) + 1)/Sampling_Freq_V;


                % First, find matches 
                TTL_Low_off = nan(size(TTL_Low, 1), 2);
                parfor z = 1:size(TTL_Low, 1)

                    clos_match = TTL_High - TTL_Low(z) < 50 ...
                        & TTL_High - TTL_Low(z) > 0;

                    if sum(clos_match) == 1
                        TTL_Low_off(z, :) = [TTL_Low(z) TTL_High(clos_match)];
                    elseif sum(clos_match) > 1
                        idx = find(clos_match);
                        if TTL_High(idx(1)) -  TTL_Low(z) < 10
                            TTL_Low_off(z, :) = [TTL_Low(z) TTL_High(idx(1))];
                        end
                    end
                end
            
                % Now ensure that there is no skipped frames
                [BinsFreq, ~] = histcounts(diff(TTL_Low_off(:, 1), 1), [0.004 0.006 0.012 0.02]);
                
                if sum(BinsFreq(2:end)) < 5
                    Camera_TTL{j} = TTL_Low_off;                    
                else
                    Camera_TTL{j} = nan;
                    disp(['Camera TTLs mismatch, check whats wrong ! ' ...
                        num2str(sum(BinsFreq(2:end))) ' TTLs are not matched'])
                end
            end
        else
            CamTTL = round(M(:, Cam));
            TTL_High = (find(diff(CamTTL,1,1) > 1) + 1)/Sampling_Freq_V;
            TTL_Low = (find(diff(CamTTL,1,1) < -1) + 1)/Sampling_Freq_V;

            try
                Camera_TTL = [TTL_High TTL_Low TTL_High-TTL_Low];
            catch
                if size(TTL_High, 1) - size(TTL_Low, 1) == 1
                    Camera_TTL = [TTL_High(1:end-1) TTL_Low];
                else
                    Camera_TTL = nan;
                    disp('Camera TTLs mismatch, check whats wrong !')
                end
            end
            
        end
    else
        Camera_TTL = nan;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Closed Loop Pump
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Pump_closed = cellfun(@(x) any(strfind(lower(x), 'pump')), Header);
    
    if any(Pump_closed)
        Subselected = M(:, Pump_closed);
        Subselected(Subselected < 0) = 0;
        
        Pump_activation_on = find((diff(Subselected) > 2) == 1);
        Pump_activation_on([false; diff(Pump_activation_on) < 10]) = [];
        
        Pump_activation_off = find((diff(Subselected) < -2) == 1); 
        Pump_activation_off([false; diff(Pump_activation_off) < 10]) = [];
        
        Pump_activity = [];
        for xx = 1:size(Pump_activation_on, 1)
            idx = find(abs((Pump_activation_off - Pump_activation_on(xx))) < 3000);
            
            if ~isempty(idx) && numel(idx) == 1
                Pump_activity = [Pump_activity; Pump_activation_on(xx) ...
                    Pump_activation_off(idx) ...
                    Pump_activation_off(idx) - Pump_activation_on(xx)];
            end
        end
        
        Pump_activity = Pump_activity/Sampling_Freq_V;
    else
        Pump_activity = nan;
    end
    
    
end

%     % Info Experiment
%     Info_experiment = struct;
%     
%     Mat_file = fullfile(Experimental_Fold, '*.mat');
%     Mat = dir(Mat_file);
%     load([Experimental_Fold '\' Mat.name]);
   %     
%     if size(Info_experiment.All_Triggers, 1) - 1 ~= size(Indices, 1)
%         error('The Experimental file doesnt match the recorded Stimulus triggers, please check!')
%     end
    
    
%     % All Triggers
%     Triggerallind = cellfun(@(x) any(strfind(x, 'All_Triggers')), Header);
%     
%     if any(All_Triggers)
%         Trigs = M(:, Triggerallind) > 0;
%         Lick_ON_OFF = [find(diff(Trigs) == 1) find(diff(Trigs) == -1) ...
%             find(diff(Trigs) == -1)-find(diff(Trigs) == 1)]/Sampling_Freq_V; % Calculate Onset, Offset and Duration
%         Lick_ON_OFF(Lick_ON_OFF(:, 3) > 0.5, :) = []; % Exclude Licks that are longer than 0.5 s
%         Licks.Lick_ON_OFF = Lick_ON_OFF;
%         Licks.Licks = Trigs;
%     else
%         Licks = [];
%     end
%     



%%% Lick Old code

        
% pd_Lick = fitdist(Frequency_Licks(Frequency_Licks < 0.3), 'Kernel', 'Width', 1/numel([0:0.005:0.3]));
% figure; plot([0:0.005:0.3], pdf(pd_Lick, [0:0.005:0.3]) * 25)
% Norm_Pois_freq = max(hist(Frequency_Licks(Frequency_Licks < 0.3), [0:0.005:0.3])) / ...
%     max(poisspdf(Frequency_Licks(Frequency_Licks < 0.3), pd_Lick.lambda));
% [f, ~] = ksdensity(Frequency_Licks(Frequency_Licks < 0.3), 'NumPoints', numel([0:0.005:0.3]));
% 
% figure; plot([0:0.005:0.3], f/numel(f)); hold on; hist(Frequency_Licks(Frequency_Licks < 0.3), [0:0.005:0.3])
% g = figure('position', [500 100 1500 1000]); 
% hist(Frequency_Licks(Frequency_Licks < 0.3), [0:0.005:0.3]); xlim([0 0.35]); hold on;
% plot(0:0.005:0.3, poisspdf(Frequency_Licks(Frequency_Licks < 0.3), pd_Lick.lambda)*Norm_Pois_freq, 'k', 'linewidth', 1.5)
% 
% 
% figure; histfit(Frequency_Licks(Frequency_Licks < 0.3), numel([0:0.005:0.3]), 'kernel')
% 
% 
% figure; hist(diff(Cl_Lick_ON_OFF(:, 1)));



