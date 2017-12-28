function [Stimulus_all_trials Stimulus_Matrix V_excite V_inhibit Cov_Matrix Cov_Matrix_raw_stimulus Stim_all_trials] = STA(p)

% This function calculates STAs, smooths the STA with cubic splining,
% calculates parameters like location of peak & trough of STA and the
% integration window of peak and integration window of trough
% In this function I call 1) compgroupSTA_estim 2) STA_parameters_significance 3) STC_excitatory_parameters_significance 4) STC_inhibitory_parameters_significance
% Created by Sudsa (20150730)
% Improved by Nima (20180101)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(p.data_dir,strcat(p.exp_id, p.dfile_suffix)));
spiketimes = eval(p.cell_id);

if isnan(p.start_TTL) || isnan(p.stop_TTL)
    A2a = A2a(1:end, 1);% this line was not the STA_simplified
else
    A2a = A2a(p.start_TTL:p.stop_TTL, 1);% this line was not the STA_simplified
end
%% ToDo: the following part of the code activated by cardinal_STA_Only_Burst should be revised
% otherwise it expected to be highly error prone
if p.cardinal_STA_Only_Burst
    % fsta spike times are recalculated so that the original spike train is burst corrected.
    % The corrected spike times can be weighted.
    % Singleton spikes can be included.
    % spiketimes(Burst_End(iter)+1) gives the first spike right after a burst.
    % Then the subsequesnt spiketimes in a burst are replaced by the first spike of the burst.
    % The number of spikes to be replaced is given by the line
    % Burst_End(iter+1)-Burst_End(iter). This is basically the length of the burst.     
    % The entire operation is being done in a for loop.
    % So in each iteration the spikes are concatenated to get the updated spiketimes
        
    diff_spiketimes = diff(spiketimes); % diff_spiketimes
    Burst_End = find(diff_spiketimes > p.ISI_gap); % p.ISI_gap is what we set based on the scatter plot of pre and post ISI
    new_spike_times = [];
 
    for iter = 1:length(Burst_End) - 1
        if (Burst_End(iter + 1) - Burst_End(iter) > 1)
            if p.weighted_burst
                new_spike_times_holder = repmat(spiketimes(Burst_End(iter) + 1), Burst_End(iter + 1) - Burst_End(iter), 1);
                new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
            elseif ~p.weighted_burst
                new_spike_times_holder = repmat(spiketimes(Burst_End(iter) + 1), 1, 1);
                new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
            end
        end
        if (Burst_End(iter + 1) - Burst_End(iter) == 1)
            if p.singleton_spikes
                new_spike_times_holder = repmat(spiketimes(Burst_End(iter) + 1), Burst_End(iter + 1) - Burst_End(iter), 1);
                new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
            end
        end
    end

    if (length(spiketimes) - Burst_End(end) > 1)
        if p.weighted_burst
            new_spike_times_holder = repmat(spiketimes(Burst_End(end) + 1), length(spiketimes) - Burst_End(end), 1); % Include all the spike times after the p.last bursting event
            new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
        end
    elseif (length(spiketimes) - Burst_End(end) == 1)
        if (~p.weighted_burst || p.singleton_spikes)
            new_spike_times_holder = spiketimes(Burst_End(end) + 1);
            new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
        end
    end
 
    if (Burst_End(1) > 1)
        if p.weighted_burst
            new_spike_times_holder = repmat(spiketimes(1), Burst_End(1), 1); % Include all the spike times before the first bursting starts
            new_spike_times = vertcat(new_spike_times_holder, new_spike_times);
        elseif ~p.weighted_burst
            new_spike_times_holder = repmat(spiketimes(1), 1, 1); % Include all the spike times before the first bursting starts
            new_spike_times = vertcat(new_spike_times_holder, new_spike_times);
        end
    elseif (Burst_End(1) == 1)
        if (p.singleton_spikes)
            new_spike_times_holder = spiketimes(1); % Include all the spike times before the first bursting starts
            new_spike_times = vertcat(new_spike_times_holder, new_spike_times);
        end
    end
    spiketimes = new_spike_times;
end

%% calculating spiketimes within each trial for a cell
Stim_all_trials = [];
Stimulus_all_trials = [];
Stimulus_Matrix = [];

NSP = 0; % Number of spikes = 0; This is the spike count variable for the STAs
inner = 0; % Total number of stimuli vectors used for STA calculation across trials

Freq = p.TTL_Pulse_Number / p.trial_length_in_secs;

p.skip_cycle = p.alternate_number * 2; % this line was not the STA_simplified

STA = zeros(2 * p.tKerLen, 1);
counter = 0;
flag_skip = true;
if isnan(p.trials_to_use)
    % Goes through trials first_trial -> last_tiral, 
    for trialIdx = p.first_trial:p.last_tiral
        if p.alternate
            if (p.alternate_number > 1) && (rem(trialIdx, p.alternate_number) == 1)
                flag_skip = xor(flag_skip,true);
            else
                flag_skip = xor(flag_skip,true);
            end
        else
            flag_skip = false;
        end
        if ~flag_skip
            counter = counter + 1;
            if p.NR
                stim_fname = strcat(fullfile(p.data_dir,'rexp_'), num2str(ceil(trialIdx)), '.txt');
                A1 = importdata(stim_fname);
            end
            Stim = [];
            % formats stimuli from text file into a vector
            for i = 8:2:p.last
                x = str2num(A1.textdata{i, 1});
                Stim = horzcat(Stim, x);
            end
            % Normalises the stimulus if specified by user
            if (p.Normalise == 1)
                Stim = (Stim - mean(Stim)) / std(Stim); %stim = Stimulus;% Sudsa Modified
            end
            ceiling = length(Stim);
            if (p.leaveit == 1)
                TTLStim = A2a((((p.TTL_Pulse_Number) * (trialIdx - 1)) + 1):((p.TTL_Pulse_Number) * (trialIdx - 1)) + p.TTL_Pulse_Number);
            elseif (p.leaveit == 2)
             
                % TTLStim = A2a((((p.TTL_Pulse_Number+p.skip)*(trialIdx-1))+(1+p.skip)):((p.TTL_Pulse_Number+p.skip)*(trialIdx)));
             
                if rem(trialIdx, p.alternate_number) == 0
                 
                    TTLStim = A2a((floor(trialIdx / p.skip_cycle) * (p.alternate_number * (p.TTL_Pulse_Number + p.skip))) + (p.TTL_Pulse_Number * (p.alternate_number - 1)) + 1:(floor(trialIdx / p.skip_cycle) * (p.alternate_number * ...
                    (p.TTL_Pulse_Number + p.skip))) + (p.TTL_Pulse_Number * p.alternate_number));
                else
                    TTLStim = A2a((floor(trialIdx / p.skip_cycle) * (p.alternate_number * (p.TTL_Pulse_Number + p.skip))) + (p.TTL_Pulse_Number * (rem(trialIdx, p.alternate_number) - 1)) + 1:(floor(trialIdx / p.skip_cycle) * (p.alternate_number * ...
                    (p.TTL_Pulse_Number + p.skip))) + (p.TTL_Pulse_Number * rem(trialIdx, p.alternate_number)));
                end
            end
         
            TTLStim = TTLStim - p.post_wait;
            tsp = spiketimes(find((spiketimes > TTLStim(1, 1) & spiketimes <= TTLStim(end, 1))), 1);
            Stim_all_trials = [Stim_all_trials Stim];
            % deletes direct RGC spikes based on lock out period
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
            % Corrects spike times based on starting TTL pulse of each trial. So all spike times will be between 0 and 100s
            tsp = tsp - TTLStim(1);
         
            TTLStim = TTLStim - TTLStim(1);
         
            TTLStim_extended = horzcat(TTLStim', TTLStim(end)+1/Freq);
            sps = histc(tsp, TTLStim_extended)';      % bin spike times
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
        counter = counter + 1;
        if p.NR
            stim_fname = strcat(fullfile(p.data_dir,'rexp_'), num2str(ceil(p.stim_to_use(counter))), '.txt');
            A1 = importdata(stim_fname);
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
        TTLStim = A2a((((p.TTL_Pulse_Number) * (I - 1)) + 1):((p.TTL_Pulse_Number) * (I - 1)) + p.TTL_Pulse_Number);
     
        TTLStim = TTLStim - p.post_wait;
        tsp = spiketimes(find((spiketimes > TTLStim(1, 1) & spiketimes <= TTLStim(end, 1))), 1);
        Stim_all_trials = [Stim_all_trials Stim];     
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
          
        tsp = tsp - TTLStim(1);
        TTLStim = TTLStim - TTLStim(1);
        TTLStim_extended = horzcat(TTLStim', TTLStim(end)+1/Freq);
        sps = histc(tsp, TTLStim_extended)';      % bin spike times
        mousetype = ' c57/bl6 wild type';
        xyz = strcat(p.cell_id, mousetype);
        title(xyz)
        [STA stimulus_matrix len NSP inner Stm_matrix] = compgroupSTA_estim(Stim, sps, p.tKerLen, STA, NSP, inner, ceiling); % compute STA
        Stimulus_all_trials = [Stimulus_all_trials Stim];
        Stimulus_Matrix = [Stimulus_Matrix; Stm_matrix];
    end
end

%% Plots begin from here
line_thickness = 2;

if p.Normalise == 1
    plt_ylim = [- 1, 1];
else
	plt_ylim = [- 1300, -300];
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
zeromarker = zeros(length(plt_ylim(1):100:plt_ylim(2)));
hold on
plot(zeromarker, plt_ylim(1):100:plt_ylim(2), 'k');
ylim([plt_ylim(1) plt_ylim(2)])
hold on
legend(['STA:spks = ', int2str(NSP)], ['Mean Stimulus : ', int2str(Freq), 'Hz stm. frq. : 35% var '])
%%
fig_names = strcat(p.names, p.cell_id);
total_time = p.trial_length_in_secs * counter;
Frequency = NSP / (counter * p.trial_length_in_secs);
xlabel({'Time (sec)'})
ylabel('Mean Stimulus (mV)')
mousetype = ' C57Bl/6 (wt) ';
figure_heading = strcat(p.year, p.cell_id, mousetype, ' FR = ', num2str(Frequency), ' Hz: ', num2str(total_time), ' sec : epiretinal', ' lock out is ', num2str(p.leave_out), 's', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes));
title(figure_heading)
set(gcf, 'name', fig_names)
fname = strcat(p.exp_id, int2str(Freq), 'Hz_', num2str(p.cont), '_', p.mean_V, '\', p.cell_id)
fname = strcat(p.work_dir, fname, '\');
if ~ exist(fname, 'dir'), mkdir(fname); end
 
    ext = '.fig';
 
    if (p.Normalise == 0)
     
        saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'not_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), ' WB is ', num2str(p.weighted_burst), ' Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
        set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
     
        saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'not_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), ' WB is ', num2str(p.weighted_burst), ' Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
     
        % close all
        figure
        %%
        y = 1 * STA;
        frame = 1 / (p.TTL_Pulse_Number / p.trial_length_in_secs);
     
        x = (- (p.tKerLen)) * (frame) + .5 * frame:frame:p.tKerLen * frame - .5 * frame;
        h = plot(x, y);
        set(gcf, 'color', 'w');
     
        saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'not_norm_for_NL', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
        set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
     
        saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'not_norm_for_NL', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
     
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
        zeromarker = zeros(length(plt_ylim(1):100:plt_ylim(2)));
        hold on
        plot(zeromarker, plt_ylim(1):100:plt_ylim(2), 'k');
        ylim([plt_ylim(1) plt_ylim(2)])
        hold on
     
        legend(['STA:spks = ', int2str(NSP)], ['Mean Stimulus : ', int2str(Freq), 'Hz stm. frq. : 35% var '])
        total_time = p.trial_length_in_secs * counter;
        Frequency = NSP / (counter * p.trial_length_in_secs);
        xlabel({'Time (sec)'})
        ylabel('Mean Stimulus (mV)')
        mousetype = ' C57Bl/6 (wt) ';
     
        figure_heading = strcat(p.year, p.cell_id, mousetype, ' FR = ', num2str(Frequency), ' Hz: ', num2str(total_time), ' sec : epiretinal', ' lock out is ', num2str(p.leave_out), 's', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes));
        title(figure_heading)
        set(gcf, 'name', fig_names)
        if ~ exist(fname, 'dir'), mkdir(fname); end
            ext = '.fig';
            saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'error_bars_not_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
            set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
            saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'error_bars_not_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
         
        else
         
            set(h, 'LineWidth', line_thickness)
            baseline = (1 * mean(Stim)) * ones(2 * p.tKerLen, 1);
            hold on
            plot(x, baseline, 'k')
            set(gcf, 'color', 'w');
            zeromarker = zeros(length(plt_ylim(1):.05:plt_ylim(2)));
            hold on
            plot(zeromarker, plt_ylim(1):.05:plt_ylim(2), 'k');
            ylim([plt_ylim(1) plt_ylim(2)])
            hold on
         
            legend(['STA:spks = ', int2str(NSP)], ['Mean Stimulus : ', int2str(Freq), 'Hz stm. frq. : 35% var '])
         
            %%
            total_time = p.trial_length_in_secs * counter;
            Frequency = NSP / (length(p.first_trial:p.last_tiral) * p.trial_length_in_secs);
            xlabel({'Time (sec)'})
            ylabel('Mean Stimulus (mV)')
            mousetype = ' C57Bl/6 (wt) ';
            set(gcf, 'color', 'w');
            saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'p.Normalised_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
            set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
         
            saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'p.Normalised_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
            figure
            y = 1 * STA;
            frame = 1 / (p.TTL_Pulse_Number / p.trial_length_in_secs);
            x = (- (p.tKerLen)) * (frame) + .5 * frame:frame:p.tKerLen * frame - .5 * frame;
            h = plot(x, y);
            set(gcf, 'color', 'w');
            saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'p.Normalised_for_NL', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
            set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
            saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'p.Normalised_for_NL', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
         
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
            plot(zeromarker, plt_ylim(1):.05:plt_ylim(2), 'k');
            errorbar(x, baseline, STA_error, 'k')
            zeromarker = zeros(length(plt_ylim(1):100:plt_ylim(2)));
            hold on
            set(gcf, 'color', 'w');
            ylim([plt_ylim(1) plt_ylim(2)])
            hold on
            legend(['STA:spks = ', int2str(NSP)], ['Mean Stimulus : ', int2str(Freq), 'Hz stm. frq. : 35% var '])
         
            %%
            total_time = p.trial_length_in_secs * counter;
            Frequency = NSP / (length(p.first_trial:p.last_tiral) * p.trial_length_in_secs);
            xlabel({'Time (sec)'})
            ylabel('Mean Stimulus (mV)')
            mousetype = ' C57Bl/6 (wt) ';
         
            figure_heading = strcat(p.year, p.cell_id, mousetype, ' FR = ', num2str(Frequency), ' Hz: ', num2str(total_time), ' sec : epiretinal', ' lock out is ', num2str(p.leave_out), 's', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes));
            title(figure_heading)
            set(gcf, 'name', fig_names)
            if ~ exist(fname, 'dir'), mkdir(fname); end
                ext = '.fig';
                saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'error_bars_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), ext], 'fig');
                set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
                saveas(gcf, [fname, p.cell_id, num2str(p.leave_out), ' ', num2str(p.first_trial), 'to', num2str(p.last_tiral), '-', 'error_bars_norm_full', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), 'WB is ', num2str(p.weighted_burst), 'Singleton is ', num2str(p.singleton_spikes), '.jpeg'], 'jpeg');
            end
            % STA parameter calculation
            STA_parameters_significance(Frequency, fname, frame, STA, Stim, p);
            % STC
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