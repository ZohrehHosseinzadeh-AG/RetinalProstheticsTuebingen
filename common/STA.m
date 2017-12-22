function [Stimulus_all_trials Stimulus_Matrix V_excite V_inhibit Cov_Matrix Cov_Matrix_raw_stimulus Stim_all_trials] = STA(p)

% This function calculates STAs, smoothes the STA with cubic splining,
% calculates parameters like location of peak & trough of STA and the
% integration window of peak and integration window of trough
% fig1 is the STA with details such as cellname, p.year of experiment, mouse
% type, firing rate, exclusion period etc
% fig2 is just the raw STA without any title or details. I use this STA for
% nonlinearity calculation.
% fig3 is the STA with errorbars
% In this function I call 1) compgroupSTA_estim
% 2) STA_parameters_significance 3) STC_excitatory_parameters_significance 4) STC_inhibitory_parameters_significance
% STA_wrapper calls this function
% Created by Sudsa (20150730)

%CODE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Stimulus_all_trials = [];

stimulus_file = fullfile(p.data_dir,'rexp_');
if p.NR == 0
    A1 = importdata(stimulus_file);
end

%% Loads mat file
inner = 0;
load(fullfile(p.data_dir,strcat(p.exp_id, p.dfile_suffix)));
cellname = eval(p.cell_id);

%% Calculates the date and name of the experiment. Later creates a folder for each experiment (based on the date)
%% in which the STAs etc are stored

%% Initialises some paramters
Stim = [];
Freq = p.TTL_Pulse_Number / p.trial_length_in_secs;
line_thickness = 2;
p.names = strcat(p.names, p.cell_id);
Stimulus_Matrix = [];
Stim_all_trials = [];
%% If STA is normalised STA, limits are reset to -1 -> 1

if p.Normalise == 1
    p.A = - 1;
    p.B = 1;
end

%% Number of spikes = 0; This is the spike count variable for the STAs
NSP = 0;
%% If you want to do fsta it recalcultes the spike times such that the original spike train is burst corrected. The corrected spike times
%% may or maynot be weighted, depending on the user input. Singleton spikes may or may not be included based on user input.
p.skip_cycle = p.alternate_number * 2; % this line was not the STA_simplified
if isnan(p.start_TTL) || isnan(p.stop_TTL)
    A2a = A2a(1:end, 1);% this line was not the STA_simplified
else
    A2a = A2a(p.start_TTL:p.stop_TTL, 1);% this line was not the STA_simplified
end

if p.cardinal_STA_Only_Burst
 
    consecutive_spike_time_gap = diff(cellname); % Calculate the spike time difference for the entire array of spike times
    Burst_End = find(consecutive_spike_time_gap > p.ISI_gap); % The p.ISI_gap is what we set based on the scatter plot of pre and post ISI. I make a list of when the bursts end in the spike train
    new_spike_times = [];
 
    for iter = 1:length(Burst_End) - 1
     
        if p.weighted_burst == 1
         
            if (Burst_End(iter + 1) - Burst_End(iter) > 1)
             
                new_spike_times_holder = repmat(cellname(Burst_End(iter) + 1), Burst_End(iter + 1) - Burst_End(iter), 1);
                new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
             
            end
         
        elseif p.weighted_burst == 0
         
            if (Burst_End(iter + 1) - Burst_End(iter) > 1)
             
                new_spike_times_holder = repmat(cellname(Burst_End(iter) + 1), 1, 1);
                new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
             
            end
         
        end
     
        if p.singleton_spikes == 1
         
            if (Burst_End(iter + 1) - Burst_End(iter) == 1)
             
                new_spike_times_holder = repmat(cellname(Burst_End(iter) + 1), Burst_End(iter + 1) - Burst_End(iter), 1);
                new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
             
            end
        end
        % cellname(Burst_End(iter)+1) gives me the first spike right after a
        % burst. I then replace the subsequesnt spike times in a burst by the
        % first spike of the burst. This is done by the above line of code. The
        % number of spikes to be replaces is given by the line
        % Burst_End(iter+1)-Burst_End(iter). This is basically the length of
        % the burst
     
        % This  entire operation is being done in a for loop. So in each iteration I concatenate the spikes
        % to get the updates spike times
     
    end
 
    if p.weighted_burst == 1
     
        if (length(cellname) - Burst_End(end) > 1)
         
            new_spike_times_holder = repmat(cellname(Burst_End(end) + 1), length(cellname) - Burst_End(end), 1); % Include all the spike times after the p.last bursting event
            new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
         
        end
     
    elseif (length(cellname) - Burst_End(end) == 1)
     
        if (p.weighted_burst == 0 || p.singleton_spikes == 1)
         
            new_spike_times_holder = cellname(Burst_End(end) + 1);
            new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
        end
     
    end
 
    if (Burst_End(1) > 1)
     
        if p.weighted_burst == 1
         
            new_spike_times_holder = repmat(cellname(1), Burst_End(1), 1); % Include all the spike times before the first bursting starts
            new_spike_times = vertcat(new_spike_times_holder, new_spike_times);
         
        elseif p.weighted_burst == 0
         
            new_spike_times_holder = repmat(cellname(1), 1, 1); % Include all the spike times before the first bursting starts
            new_spike_times = vertcat(new_spike_times_holder, new_spike_times);
         
        end
     
    elseif (Burst_End(1) == 1)
     
        if (p.singleton_spikes == 1)
         
            new_spike_times_holder = cellname(1); % Include all the spike times before the first bursting starts
            new_spike_times = vertcat(new_spike_times_holder, new_spike_times);
         
        end
     
    end
 
    cellname = new_spike_times;
 
end

%%

STA = zeros(2 * p.tKerLen, 1);
abc = 0;

counter = 0;

if isnan(p.trials_to_use)
 
    %% Goes through trials ii -> jj, calculating the spike times in each trial for a cell.
    for I = p.ii:p.jj
     
        %%not present in the simplified STA
        if p.alternate
         
            if I > 1
             
                if p.alternate_number > 1
                 
                    if (rem(I, p.alternate_number) == 1)
                     
                        if p.flag_skip
                         
                            p.flag_skip = 0;
                         
                        elseif p.flag_skip == 0
                         
                            p.flag_skip = 1;
                         
                        end
                     
                    end
                 
                elseif p.alternate_number < 1
                 
                    if p.flag_skip
                     
                        p.flag_skip = 0;
                     
                    elseif p.flag_skip == 0
                     
                        p.flag_skip = 1;
                     
                    end
                 
                end
             
            end
            %%
        elseif p.alternate == 0
         
            p.flag_skip = 1;
        end
     
        if p.flag_skip
         
            I
            counter = counter + 1;
            if p.NR == 0
                p.ii = 1;
            else %% loads stimuli for each trial, if the experiment has non repeating trials
                p.ii = I;
                name = strcat(stimulus_file, num2str(ceil(p.ii)), '.txt');
                name
                A1 = importdata(name);
            end
         
            Stim = [];
            %% formats stimuli from text file into a vector
            for i = 8:2:p.last
                x = str2num(A1.textdata{i, 1});
                Stim = horzcat(Stim, x);
            end
         
            %% Normalises the stimulus if specified by user
            if (p.Normalise == 1)
                Stim = (Stim - mean(Stim)) / std(Stim); %stim = Stimulus;% Sudsa Modified
            end
            ceiling = length(Stim);
            if (p.leaveit == 1)
                I
                TTLStim = A2a((((p.TTL_Pulse_Number) * (I - 1)) + 1):((p.TTL_Pulse_Number) * (I - 1)) + p.TTL_Pulse_Number);
             
            elseif (p.leaveit == 2)
             
                %             TTLStim = A2a((((p.TTL_Pulse_Number+p.skip)*(I-1))+(1+p.skip)):((p.TTL_Pulse_Number+p.skip)*(I)));
             
                if rem(I, p.alternate_number) == 0
                 
                    TTLStim = A2a((floor(I / p.skip_cycle) * (p.alternate_number * (p.TTL_Pulse_Number + p.skip))) + (p.TTL_Pulse_Number * (p.alternate_number - 1)) + 1:(floor(I / p.skip_cycle) * (p.alternate_number * ...
                    (p.TTL_Pulse_Number + p.skip))) + (p.TTL_Pulse_Number * p.alternate_number));
                 
                else
                 
                    TTLStim = A2a((floor(I / p.skip_cycle) * (p.alternate_number * (p.TTL_Pulse_Number + p.skip))) + (p.TTL_Pulse_Number * (rem(I, p.alternate_number) - 1)) + 1:(floor(I / p.skip_cycle) * (p.alternate_number * ...
                    (p.TTL_Pulse_Number + p.skip))) + (p.TTL_Pulse_Number * rem(I, p.alternate_number)));
                 
                end
             
            end
         
            TTLStim = TTLStim - p.post_wait;
            tsp = cellname(find((cellname > TTLStim(1, 1) & cellname <= TTLStim(end, 1))), 1);
            Stim_all_trials = [Stim_all_trials Stim];
            %%
            %% deletes direct RGC spikes based on lock out period
         
            if p.leave_out > 0
                tsp_del_loc = [];
                for count = 1:length(TTLStim)
                 
                    tsp_del_temp_loc = find(tsp >= TTLStim(count) & tsp <= TTLStim(count) + p.leave_out);
                    [m n] = size(tsp_del_temp_loc);
                 
                    if (m > 1)
                        tsp_del_temp_loc = tsp_del_temp_loc';
                    end
                    tsp_del_loc = horzcat(tsp_del_loc, tsp_del_temp_loc);
                end
             
                tsp_before = tsp;
                tsp(tsp_del_loc) = [];
            end
         
            %%    Corrects spike times based on starting TTL pulse of each trial. So all spike times will be between 0 and 100s
         
            tsp = tsp - TTLStim(1);
         
            TTLStim = TTLStim - TTLStim(1);
         
            TTLStim_extended = horzcat(TTLStim', TTLStim(end)+1/Freq);
            sps = histc(tsp, TTLStim_extended)';      % bin spike times
            abc = abc + sum(sps);
            mousetype = ' c57/bl6 wild type';
            xyz = strcat(p.cell_id, mousetype);
            title(xyz)
         
            [STA stimulus_matrix len NSP inner Stm_matrix] = compgroupSTA_estim(Stim, sps, p.tKerLen, STA, NSP, inner, ceiling); % compute STA
         
            Stimulus_all_trials = [Stimulus_all_trials Stim];
         
            Stimulus_Matrix = [Stimulus_Matrix; Stm_matrix];
         
        end
    end
 
else
 
    for J = 1:length(p.trials_to_use)
     
        I = p.trials_to_use(J)
     
        %%
     
        counter = counter + 1;
     
        if p.NR == 0
            p.ii = 1;
        else
            p.ii = I;
            name = strcat(stimulus_file, num2str(ceil(p.stim_to_use(counter))), '.txt');
            A1 = importdata(name);
        end
     
        Stim = [];
        for i = 8:2:p.last
            x = str2num(A1.textdata{i, 1});
            Stim = horzcat(Stim, x);
        end
     
        if (p.Normalise == 1)
            Stim = (Stim - mean(Stim)) / std(Stim); %stim = Stimulus;% Sudsa Modified
        end
        ceiling = length(Stim);
     
        size(A2a)
     
        TTLStim = A2a((((p.TTL_Pulse_Number) * (I - 1)) + 1):((p.TTL_Pulse_Number) * (I - 1)) + p.TTL_Pulse_Number);
     
        TTLStim = TTLStim - p.post_wait;
        tsp = cellname(find((cellname > TTLStim(1, 1) & cellname <= TTLStim(end, 1))), 1);
        Stim_all_trials = [Stim_all_trials Stim];
        %%
     
        if p.leave_out > 0
            tsp_del_loc = [];
            for count = 1:length(TTLStim)
             
                tsp_del_temp_loc = find(tsp >= TTLStim(count) & tsp <= TTLStim(count) + p.leave_out);
                [m n] = size(tsp_del_temp_loc);
             
                if (m > 1)
                    tsp_del_temp_loc = tsp_del_temp_loc';
                end
                tsp_del_loc = horzcat(tsp_del_loc, tsp_del_temp_loc);
            end
         
            tsp_before = tsp;
            tsp(tsp_del_loc) = [];
        end
     
        %%
     
        tsp = tsp - TTLStim(1);
     
        TTLStim = TTLStim - TTLStim(1);
     
        TTLStim_extended = horzcat(TTLStim', TTLStim(end)+1/Freq);
        sps = histc(tsp, TTLStim_extended)';      % bin spike times
        abc = abc + sum(sps);
        mousetype = ' c57/bl6 wild type';
        xyz = strcat(p.cell_id, mousetype);
        title(xyz)
     
        [STA stimulus_matrix len NSP inner Stm_matrix] = compgroupSTA_estim(Stim, sps, p.tKerLen, STA, NSP, inner, ceiling); % compute STA
     
        Stimulus_all_trials = [Stimulus_all_trials Stim];
     
        Stimulus_Matrix = [Stimulus_Matrix; Stm_matrix];
     
    end
 
end

STA = STA / NSP; % Divides the Stimulus Sum/ Total Number of Spikes

y = 1 * STA;

frame = 1 / (p.TTL_Pulse_Number / p.trial_length_in_secs);

x = (- (p.tKerLen)) * (frame) + .5 * frame:frame:p.tKerLen * frame - .5 * frame;
h = plot(x, y);
set(h, 'LineWidth', line_thickness)
baseline = (1 * mean(Stim)) * ones(2 * p.tKerLen, 1);
hold on
plot(x, baseline, 'k')
set(gcf, 'color', 'w');
zeromarker = zeros(length(p.A:100:p.B));
hold on
plot(zeromarker, p.A:100:p.B, 'k');
ylim([p.A p.B])
hold on

legend(['STA:spks = ', int2str(NSP)], ['Mean Stimulus : ', int2str(Freq), 'Hz stm. frq. : 35% var '])

%%

total_time = p.trial_length_in_secs * counter;
Frequency = NSP / (counter * p.trial_length_in_secs);
xlabel({'Time (sec)'})
ylabel('Mean Stimulus (mV)')
mousetype = ' C57Bl/6 (wt) ';

figure_heading = strcat(p.year, p.cell_id, mousetype, ' FR = ', num2str(Frequency), ' Hz: ', num2str(total_time), ' sec : epiretinal', ' lock out is ', num2str(p.leave_out), 's', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes));
title(figure_heading)

set(gcf, 'name', p.names)

fname = strcat(p.exp_id, int2str(Freq), 'Hz_', num2str(p.cont), '_', p.mean_V, '\', p.cell_id)

fname = strcat(p.work_dir, fname, '\');

dump = fname

if ~ exist(dump, 'dir'), mkdir(dump); end
 
    ext = '.fig';
 
    if (p.Normalise == 0)
     
        saveas(gcf, [dump, p.cell_id, num2str(p.leave_out), ' ', num2str(p.ii), 'to', num2str(p.jj), '-', 'not_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), ' WB is ', num2str(p.weighted_burst), ' Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
        set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
     
        saveas(gcf, [dump, p.cell_id, num2str(p.leave_out), ' ', num2str(p.ii), 'to', num2str(p.jj), '-', 'not_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), ' WB is ', num2str(p.weighted_burst), ' Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
     
        % close all
        figure
        %%
        y = 1 * STA;
        frame = 1 / (p.TTL_Pulse_Number / p.trial_length_in_secs);
     
        x = (- (p.tKerLen)) * (frame) + .5 * frame:frame:p.tKerLen * frame - .5 * frame;
        h = plot(x, y);
        set(gcf, 'color', 'w');
     
        saveas(gcf, [dump, p.cell_id, num2str(p.leave_out), ' ', num2str(p.ii), 'to', num2str(p.jj), '-', 'not_norm_for_NL', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
        set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
     
        saveas(gcf, [dump, p.cell_id, num2str(p.leave_out), ' ', num2str(p.ii), 'to', num2str(p.jj), '-', 'not_norm_for_NL', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
     
        %% error bars
     
        figure
        y = 1 * STA;
        STA_error = std(y(length(y) / 2 + 1:end)) / sqrt(p.tKerLen);
        STA_error = STA_error * ones(length(y), 1);
        dim = max(size(Stimulus_Matrix));
        for ijk = 1:2 * p.tKerLen
            errorbar_sta(ijk) = std(Stimulus_Matrix(:, ijk)) / sqrt(dim);
        end
     
        frame = 1 / (p.TTL_Pulse_Number / p.trial_length_in_secs);
        x = (- (p.tKerLen)) * (frame) + .5 * frame:frame:p.tKerLen * frame - .5 * frame;
        h = errorbar(x, y, errorbar_sta);
        set(h, 'LineWidth', line_thickness)
        baseline = (1 * mean(Stim)) * ones(2 * p.tKerLen, 1);
        hold on
     
        errorbar(x, baseline, STA_error, 'k')
        set(gcf, 'color', 'w');
        zeromarker = zeros(length(p.A:100:p.B));
        hold on
        plot(zeromarker, p.A:100:p.B, 'k');
        ylim([p.A p.B])
        hold on
     
        legend(['STA:spks = ', int2str(NSP)], ['Mean Stimulus : ', int2str(Freq), 'Hz stm. frq. : 35% var '])
        total_time = p.trial_length_in_secs * counter;
        Frequency = NSP / (counter * p.trial_length_in_secs);
        xlabel({'Time (sec)'})
        ylabel('Mean Stimulus (mV)')
        mousetype = ' C57Bl/6 (wt) ';
     
        figure_heading = strcat(p.year, p.cell_id, mousetype, ' FR = ', num2str(Frequency), ' Hz: ', num2str(total_time), ' sec : epiretinal', ' lock out is ', num2str(p.leave_out), 's', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes));
        title(figure_heading)
        set(gcf, 'name', p.names)
        if ~ exist(dump, 'dir'), mkdir(dump); end
            ext = '.fig';
            saveas(gcf, [dump, p.cell_id, num2str(p.leave_out), ' ', num2str(p.ii), 'to', num2str(p.jj), '-', 'error_bars_not_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
            set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
            saveas(gcf, [dump, p.cell_id, num2str(p.leave_out), ' ', num2str(p.ii), 'to', num2str(p.jj), '-', 'error_bars_not_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
         
        else
         
            set(h, 'LineWidth', line_thickness)
            baseline = (1 * mean(Stim)) * ones(2 * p.tKerLen, 1);
            hold on
            plot(x, baseline, 'k')
            set(gcf, 'color', 'w');
            zeromarker = zeros(length(p.A:.05:p.B));
            hold on
            plot(zeromarker, p.A:.05:p.B, 'k');
            ylim([p.A p.B])
            hold on
         
            legend(['STA:spks = ', int2str(NSP)], ['Mean Stimulus : ', int2str(Freq), 'Hz stm. frq. : 35% var '])
         
            %%
            total_time = p.trial_length_in_secs * counter;
            Frequency = NSP / (length(p.ii:p.jj) * p.trial_length_in_secs);
            xlabel({'Time (sec)'})
            ylabel('Mean Stimulus (mV)')
            mousetype = ' C57Bl/6 (wt) ';
            set(gcf, 'color', 'w');
            saveas(gcf, [dump, p.cell_id, num2str(p.leave_out), ' ', num2str(p.ii), 'to', num2str(p.jj), '-', 'p.Normalised_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
            set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
         
            saveas(gcf, [dump, p.cell_id, num2str(p.leave_out), ' ', num2str(p.ii), 'to', num2str(p.jj), '-', 'p.Normalised_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
            figure
            y = 1 * STA;
            frame = 1 / (p.TTL_Pulse_Number / p.trial_length_in_secs);
            x = (- (p.tKerLen)) * (frame) + .5 * frame:frame:p.tKerLen * frame - .5 * frame;
            h = plot(x, y);
            set(gcf, 'color', 'w');
            saveas(gcf, [dump, p.cell_id, num2str(p.leave_out), ' ', num2str(p.ii), 'to', num2str(p.jj), '-', 'p.Normalised_for_NL', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
            set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
            saveas(gcf, [dump, p.cell_id, num2str(p.leave_out), ' ', num2str(p.ii), 'to', num2str(p.jj), '-', 'p.Normalised_for_NL', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
         
            %% error bars
            y = 1 * STA;
            STA_error = std(y(length(y) / 2 + 1:end)) / sqrt(p.tKerLen);
            STA_error = STA_error * ones(length(y), 1);
            frame = 1 / (p.TTL_Pulse_Number / p.trial_length_in_secs);
            x = (- (p.tKerLen)) * (frame) + .5 * frame:frame:p.tKerLen * frame - .5 * frame;
            h = plot(x, y);
            set(h, 'LineWidth', line_thickness)
            baseline = (1 * mean(Stim)) * ones(2 * p.tKerLen, 1);
            hold on
            plot(zeromarker, p.A:.05:p.B, 'k');
            errorbar(x, baseline, STA_error, 'k')
            zeromarker = zeros(length(p.A:100:p.B));
            hold on
            set(gcf, 'color', 'w');
            ylim([p.A p.B])
            hold on
            legend(['STA:spks = ', int2str(NSP)], ['Mean Stimulus : ', int2str(Freq), 'Hz stm. frq. : 35% var '])
         
            %%
            total_time = p.trial_length_in_secs * counter;
            Frequency = NSP / (length(p.ii:p.jj) * p.trial_length_in_secs);
            xlabel({'Time (sec)'})
            ylabel('Mean Stimulus (mV)')
            mousetype = ' C57Bl/6 (wt) ';
         
            figure_heading = strcat(p.year, p.cell_id, mousetype, ' FR = ', num2str(Frequency), ' Hz: ', num2str(total_time), ' sec : epiretinal', ' lock out is ', num2str(p.leave_out), 's', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes));
            title(figure_heading)
            set(gcf, 'name', p.names)
            if ~ exist(dump, 'dir'), mkdir(dump); end
                ext = '.fig';
                saveas(gcf, [dump, p.cell_id, num2str(p.leave_out), ' ', num2str(p.ii), 'to', num2str(p.jj), '-', 'error_bars_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
                set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
                saveas(gcf, [dump, p.cell_id, num2str(p.leave_out), ' ', num2str(p.ii), 'to', num2str(p.jj), '-', 'error_bars_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
             
            end
            %% STA parameter calculation
            STA_parameters_significance(Frequency, dump, frame, STA, Stim, p);
            %% STC
         
            if p.STC_Analysis
                [V_inhibit V_excite Cov_Matrix Cov_Matrix_raw_stimulus mu sd] = STC_excitatory_parameters_significance(p.tKerLen, frame, STA, Stim_all_trials, Stimulus_Matrix, p.cell_id, p.year, p.cardinal_STA_Only_Burst);
                %  [V_inhibit Cov_Matrix Cov_Matrix_raw_stimulus Stim_all_trials] = STC_inhibitory_parameters_significance(p.tKerLen, frame, STA, Stim_all_trials, Stimulus_Matrix, p.cell_id, p.year, p.cardinal_STA_Only_Burst, mu, sd);
             
            else
             
                V_excite = 0;
                V_inhibit = 0;
                Cov_Matrix = 0;
                Cov_Matrix_raw_stimulus = 0;
             
            end
         
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     