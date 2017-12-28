% This function calculates calculates parameters like location of peak & trough of STA and the
% integration window of peak and integration window of trough
% fig4 is the cubic splined STA. Marked with black circles, are the width of
% the integration window of peak and width of integration windown of trough
% STA calls this function
% Created by Sudsa (20150730)

%% paramters of the function call
function STA_parameters_significance(Frequency, dump, frame, STA, Stim, p)
% Frequency is the frequency of stimulation
% dump is the location to which the figures should be saved. Its calculated
% implicitly in the previous code STA.  p.tKerLen is 25. frame is .04s.
% STA is the STA vector.
% Stim is the stimulus of the last stimulus trial. I use it to calculate baseline mean and variance. I use these values for significance calculation later.
% p.cell_id is the name of the cell. le
% p.leave_out is the exclusion period of the STA.
% p.ii and p.jj are the trial number start and end.
% p.year is the p.year of the experiment
% p.cardinal_STA_Only_Burst is 1 if you do fSTA
% p.vstim is 0 if this function is called by estim script. Else it is 1 if it is called by p.vstim script
% p.weighted_burst. It is 1 if you do weighted burst
% p.singleton_spikes. It is 1 if you include singleton spikes

%% CODE
if p.Normalise == 1
    plt_ylim = [- 1, 1];
else
	plt_ylim = [- 1300, -300];
end

%% cubic splined sta
t = (- (p.tKerLen)) * (frame) + .5 * frame:frame:p.tKerLen * frame - .5 * frame;

tt = STA;
ttt = t(1):.001:t(end);
if p.single_pulse_activation_correction
    
    tt(length(STA) / 2) = mean(Stim);
    
end
T = spline(t, tt, ttt);
figure
plot(ttt, T, 'LineWidth', 2)
ylim([plt_ylim(1) plt_ylim(2)])
set(gcf, 'color', 'w');
baseline = (1 * mean(Stim)) * ones(length(T), 1);

hold on
plot(ttt, baseline, 'k')
step_size = (plt_ylim(2) - plt_ylim(1)) / 10;
zeromarker = zeros(length(plt_ylim(1):step_size:plt_ylim(2)));
plot(zeromarker, plt_ylim(1):step_size:plt_ylim(2), 'k');

%% significance of peak and trough of STA

baseline_sta = mean(STA(length(STA) / 2 + 1:end));
std_baselines = std(STA(length(STA) / 2 + 1:end));

alphaville = 1 - nthroot(.95, p.tKerLen);

alphaville = alphaville / 2;

peak_sta = min(STA(1:length(STA) / 2));
[H_peak P_peak] = ztest(peak_sta, baseline_sta, std_baselines, alphaville);

trough_sta = max(STA(1:length(STA) / 2));
[H_trough P_trough] = ztest(trough_sta, baseline_sta, std_baselines, alphaville);

%% Location of peak and trough and their values

if ~ p.vstim
    
    peak_T = min(T(1:(length(T) / 2)));
    
elseif p.vstim
    
    peak_T = min(T(1:(length(T) / 2)));
    
end

% trough_T = max(T(length(T)/4:length(T)/2));
trough_T = max(T(1:length(T) / 2));

baseline_avg_right = mean(T(1 + (length(T) / 2):length(T)));
baseline_std_right = std(T(1 + (length(T) / 2):length(T)));

trough_location = find(T(1:length(T) / 2) == trough_T);
peak_location = find(T(1:length(T) / 2) == peak_T);

if ~ p.vstim
    trough_location_time = abs(ttt(1) + trough_location * .001 - .001);
    peak_location_time = abs(ttt(1) + peak_location * .001 - .001);
    
elseif p.vstim
    
    trough_location_time = abs(ttt(1) + peak_location * .001 - .001);
    peak_location_time = abs(ttt(1) + trough_location * .001 - .001);
    
end

%% integration time calculation peak

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if H_peak == 1
    
    locs = find(T(1:length(T) / 2) < mean(Stim));
    i_zerocrossing = (locs(find(locs == peak_location)));
    j_zerocrossing = (locs(find(locs == peak_location)));
    i_peak = i_zerocrossing;
    
else
    
    i_peak = - 1000;
    
end

if H_trough == 1
    
    locs = find(T(1:length(T) / 2) > mean(Stim));
    i_zerocrossing = (locs(find(locs == trough_location)));
    j_zerocrossing = (locs(find(locs == trough_location)));
    i_trough = i_zerocrossing;
    
else
    
    i_trough = - 1000;
    
end

if (i_peak > i_trough)
    
    locs = find(T(1:length(T) / 2) < mean(Stim));
    i_zerocrossing = (locs(find(locs == peak_location)));
    j_zerocrossing = (locs(find(locs == peak_location)));
    
    i_sig = i_zerocrossing;
    j_sig = j_zerocrossing;
    
    while (T(i_zerocrossing) < mean(Stim))
        
        i_zerocrossing = i_zerocrossing - 1;
        if (i_zerocrossing == 0)
            i_zerocrossing = 1;
            break;
        end
        
    end
    
    while ((T(i_sig) < mean(Stim)) && ztest(T(i_sig), baseline_sta, std_baselines, alphaville))
        
        i_sig = i_sig - 1;
        if (i_sig == 0)
            i_sig = 1;
            break;
        end
        
    end
    
    plot(t(1) - .001 + (i_zerocrossing * .001), T(i_zerocrossing), 'ok', 'LineWidth', 2)
    hold on
    plot(t(1) - .001 + (i_sig * .001), T(i_sig), 'or', 'LineWidth', 2)
    
    while (T(j_zerocrossing) < mean(Stim))
        j_zerocrossing = j_zerocrossing + 1;
        if (j_zerocrossing > length(T) / 2)
            break;
        end
        
    end
    
    while ((T(j_sig) < mean(Stim)) && ztest(T(j_sig), baseline_sta, std_baselines, alphaville))
        
        j_sig = j_sig + 1;
        if (j_sig > length(T) / 2)
            break;
        end
        
    end
    
    plot(t(1) - .001 + (j_zerocrossing * .001), T(j_zerocrossing), 'ok', 'LineWidth', 2)
    hold on
    plot(t(1) - .001 + (j_sig * .001), T(j_sig), 'or', 'LineWidth', 2)
    
    if ~ p.vstim
        
        peak_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
        peak_Integration_time_sig = (j_sig - i_sig) * .001;
        
    elseif p.vstim
        
        trough_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
        trough_Integration_time_sig = (j_sig - i_sig) * .001;
        
    end
    
    i_zerocrossing_rebound = i_zerocrossing;
    j_zerocrossing_rebound_end = i_zerocrossing_rebound;
    
    while (T(i_zerocrossing_rebound) > mean(Stim))
        i_zerocrossing_rebound = i_zerocrossing_rebound - 1;
        if (i_zerocrossing_rebound == 1)
            break;
        end
    end
    
    plot(t(1) - .001 + (i_zerocrossing_rebound * .001), T(i_zerocrossing_rebound), 'ok', 'LineWidth', 2)
    plot(t(1) - .001 + ((i_zerocrossing + 1) * .001), T(i_zerocrossing + 1), 'ok', 'LineWidth', 2)
    
    trough_amp_rebound = max(T(i_zerocrossing_rebound:j_zerocrossing_rebound_end));
    
    [H_trough_rebound P_trough_rebound] = ztest(trough_amp_rebound, baseline_sta, std_baselines, alphaville);
    
    if H_trough_rebound
        
        if ~ p.vstim
            trough_Integration_time = (j_zerocrossing_rebound_end - i_zerocrossing_rebound) * .001;
            
        elseif p.vstim
            
            peak_Integration_time = (j_zerocrossing_rebound_end - i_zerocrossing_rebound) * .001;
            
        end
        
        trough_location_rebound = find(T(1:length(T) / 2) == trough_amp_rebound);
        trough_location_time_rebound = abs(ttt(1) + trough_location_rebound * .001 - .001);
        peak_location_time_rebound = - 10;
        
    else
        
        if ~ p.vstim
            
            trough_Integration_time = 0;
            peak_location_time_rebound = - 10;
            trough_location_time_rebound = - 10;
        elseif p.vstim
            
            peak_Integration_time = 0;
            
        end
        
    end
    
    if H_trough == 1
        
        locs = find(T(1:length(T) / 2) > mean(Stim));
        i_zerocrossing = (locs(find(locs == trough_location)));
        j_zerocrossing = (locs(find(locs == trough_location)));
        i_sig = i_zerocrossing;
        j_sig = j_zerocrossing;
        
        while ((T(i_sig) > mean(Stim)) && ztest(T(i_sig), baseline_sta, std_baselines, alphaville))
            
            i_sig = i_sig - 1;
            if (i_sig == 0)
                
                i_sig = 1;
                break;
            end
            
        end
        
        while ((T(j_sig) > mean(Stim)) && ztest(T(j_sig), baseline_sta, std_baselines, alphaville))
            
            j_sig = j_sig + 1;
            if (j_sig > length(T) / 2)
                break;
            end
            
        end
        
        trough_Integration_time_sig = (j_sig - i_sig) * .001;
        
        plot(t(1) - .001 + (j_sig * .001), T(j_sig), 'or', 'LineWidth', 2)
        hold on
        plot(t(1) - .001 + (i_sig * .001), T(i_sig), 'or', 'LineWidth', 2)
        
    else
        
        trough_Integration_time_sig = 0;
        
    end
    
elseif (i_peak < i_trough)
    
    locs = find(T(1:length(T) / 2) > mean(Stim));
    i_zerocrossing = (locs(find(locs == trough_location)));
    j_zerocrossing = (locs(find(locs == trough_location)));
    i_sig = i_zerocrossing;
    j_sig = j_zerocrossing;
    
    while (T(i_zerocrossing) > mean(Stim))
        
        i_zerocrossing = i_zerocrossing - 1;
        
        if (i_zerocrossing == 0)
            i_zerocrossing = 1;
            break;
            
        end
        
    end
    
    while ((T(i_sig) > mean(Stim)) && ztest(T(i_sig), baseline_sta, std_baselines, alphaville))
        
        i_sig = i_sig - 1;
        if (i_sig == 0)
            i_sig = 1;
            break;
        end
        
    end
    
    plot(t(1) - .001 + (i_zerocrossing * .001), T(i_zerocrossing), 'ok', 'LineWidth', 2)
    hold on
    plot(t(1) - .001 + (i_sig * .001), T(i_sig), 'or', 'LineWidth', 2)
    
    while (T(j_zerocrossing) > mean(Stim))
        j_zerocrossing = j_zerocrossing + 1;
        if (j_zerocrossing > length(T) / 2)
            break;
        end
    end
    
    while ((T(j_sig) > mean(Stim)) && ztest(T(j_sig), baseline_sta, std_baselines, alphaville))
        
        j_sig = j_sig + 1;
        if (j_sig > length(T) / 2)
            break;
        end
    end
    
    plot(t(1) - .001 + (j_zerocrossing * .001), T(j_zerocrossing), 'ok', 'LineWidth', 2)
    hold on
    plot(t(1) - .001 + (j_sig * .001), T(j_sig), 'or', 'LineWidth', 2)
    
    if ~ p.vstim
        
        trough_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
        trough_Integration_time_sig = (j_sig - i_sig) * .001;
        
    elseif p.vstim
        
        peak_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
        peak_Integration_time_sig = (j_sig - i_sig) * .001;
        
    end
    
    i_zerocrossing_rebound = i_zerocrossing;
    j_zerocrossing_rebound_end = i_zerocrossing_rebound;
    
    while (T(i_zerocrossing_rebound) < mean(Stim))
        i_zerocrossing_rebound = i_zerocrossing_rebound - 1;
        if (i_zerocrossing_rebound == 1)
            break;
        end
    end
    
    plot(t(1) - .001 + (i_zerocrossing_rebound * .001), T(i_zerocrossing_rebound), 'ok', 'LineWidth', 2)
    plot(t(1) - .001 + ((i_zerocrossing + 1) * .001), T(i_zerocrossing + 1), 'ok', 'LineWidth', 2)
    
    peak_amp_rebound = min(T(i_zerocrossing_rebound:j_zerocrossing_rebound_end));
    
    [H_peak_rebound P_peak_rebound] = ztest(peak_amp_rebound, baseline_sta, std_baselines, alphaville);
    
    if H_peak_rebound
        
        if ~ p.vstim
            peak_Integration_time = (j_zerocrossing_rebound_end - i_zerocrossing_rebound) * .001;
            
        elseif p.vstim
            
            trough_Integration_time = (j_zerocrossing_rebound_end - i_zerocrossing_rebound) * .001;
            
        end
        
        peak_location_rebound = find(T(1:length(T) / 2) == peak_amp_rebound);
        peak_location_time_rebound = abs(ttt(1) + peak_location_rebound * .001 - .001);
        trough_location_time_rebound = - 10;
        
    else
        
        if ~ p.vstim
            
            peak_Integration_time = 0;
            peak_location_time_rebound = - 10;
            trough_location_time_rebound = - 10;
            
        elseif p.vstim
            
            trough_Integration_time = 0;
            
        end
        
    end
    
    if H_peak == 1
        
        locs = find(T(1:length(T) / 2) < mean(Stim));
        i_zerocrossing = (locs(find(locs == peak_location)));
        j_zerocrossing = (locs(find(locs == peak_location)));
        i_sig = i_zerocrossing;
        j_sig = j_zerocrossing;
        
        while ((T(i_sig) < mean(Stim)) && ztest(T(i_sig), baseline_sta, std_baselines, alphaville))
            
            i_sig = i_sig - 1;
            if (i_sig == 0)
                
                i_sig = 1;
                break;
            end
            
        end
        
        while ((T(j_sig) < mean(Stim)) && ztest(T(j_sig), baseline_sta, std_baselines, alphaville))
            
            j_sig = j_sig + 1;
            if (j_sig > length(T) / 2)
                break;
            end
            
        end
        
        peak_Integration_time_sig = (j_sig - i_sig) * .001;
        plot(t(1) - .001 + (j_sig * .001), T(j_sig), 'or', 'LineWidth', 2)
        hold on
        plot(t(1) - .001 + (i_sig * .001), T(i_sig), 'or', 'LineWidth', 2)
        
    else
        
        peak_Integration_time_sig = 0;
        
    end
    
elseif (i_peak == i_trough)
    
    peak_Integration_time = 0;
    trough_Integration_time = 0;
    
end

%% p.single_pulse_activation_correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~ p.vstim
    
    if p.single_pulse_activation_correction
        
        if ((((peak_location_time < .04) && (H_peak == 1)) || (((trough_location_time < .04) && (H_trough == 1)))))
            disp('single pulse')
            
            tt(25) = mean(baseline);
            figure;
            title('check')
            T = spline(t, tt, ttt);
            
            plot(ttt, T, 'LineWidth', 2)
            ylim([plt_ylim(1) plt_ylim(2)])
            set(gcf, 'color', 'w');
            baseline = (1 * mean(Stim)) * ones(length(T), 1);
            hold on
            plot(ttt, baseline, 'k')
            step_size = (plt_ylim(2) - plt_ylim(1)) / 10;
            zeromarker = zeros(length(plt_ylim(1):step_size:plt_ylim(2)));
            plot(zeromarker, plt_ylim(1):step_size:plt_ylim(2), 'k');
            
            %% significance of peak and trough of STA
            
            baseline_sta = mean(STA(length(STA) / 2 + 1:end));
            std_baselines = std(STA(length(STA) / 2 + 1:end));
            
            alphaville = 1 - nthroot(.95, p.tKerLen);
            
            alphaville = alphaville / 2;
            
            peak_sta = min(STA(1:length(STA) / 2));
            [H_peak P_peak] = ztest(peak_sta, baseline_sta, std_baselines, alphaville);
            
            trough_sta = max(STA(1:length(STA) / 2));
            [H_trough P_trough] = ztest(trough_sta, baseline_sta, std_baselines, alphaville);
            
            %% Location of peak and trough and their values
            
            if ~ p.vstim
                
                peak_T = min(T(1:(length(T) / 2)));
                
            elseif p.vstim
                
                peak_T = min(T(1:(length(T) / 2)));
                
            end
            
            % trough_T = max(T(length(T)/4:length(T)/2));
            trough_T = max(T(1:length(T) / 2));
            
            baseline_avg_right = mean(T(1 + (length(T) / 2):length(T)));
            baseline_std_right = std(T(1 + (length(T) / 2):length(T)));
            
            trough_location = find(T(1:length(T) / 2) == trough_T);
            peak_location = find(T(1:length(T) / 2) == peak_T);
            
            if ~ p.vstim
                trough_location_time = abs(ttt(1) + trough_location * .001 - .001);
                peak_location_time = abs(ttt(1) + peak_location * .001 - .001);
                
            elseif p.vstim
                
                trough_location_time = abs(ttt(1) + peak_location * .001 - .001);
                peak_location_time = abs(ttt(1) + trough_location * .001 - .001);
                
            end
            
            %% integration time calculation peak
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            if H_peak == 1
                
                locs = find(T(1:length(T) / 2) < mean(Stim));
                i_zerocrossing = (locs(find(locs == peak_location)));
                j_zerocrossing = (locs(find(locs == peak_location)));
                i_peak = i_zerocrossing;
                
            else
                
                i_peak = - 1000;
                
            end
            
            if H_trough == 1
                
                locs = find(T(1:length(T) / 2) > mean(Stim));
                i_zerocrossing = (locs(find(locs == trough_location)));
                j_zerocrossing = (locs(find(locs == trough_location)));
                i_trough = i_zerocrossing;
                
            else
                
                i_trough = - 1000;
                
            end
            
            if (i_peak > i_trough)
                
                locs = find(T(1:length(T) / 2) < mean(Stim));
                i_zerocrossing = (locs(find(locs == peak_location)));
                j_zerocrossing = (locs(find(locs == peak_location)));
                
                i_sig = i_zerocrossing;
                j_sig = j_zerocrossing;
                
                while (T(i_zerocrossing) < mean(Stim))
                    
                    i_zerocrossing = i_zerocrossing - 1;
                    if (i_zerocrossing == 0)
                        i_zerocrossing = 1;
                        break;
                    end
                    
                end
                
                while ((T(i_sig) < mean(Stim)) && ztest(T(i_sig), baseline_sta, std_baselines, alphaville))
                    
                    i_sig = i_sig - 1;
                    if (i_sig == 0)
                        i_sig = 1;
                        break;
                    end
                    
                end
                
                plot(t(1) - .001 + (i_zerocrossing * .001), T(i_zerocrossing), 'ok', 'LineWidth', 2)
                hold on
                plot(t(1) - .001 + (i_sig * .001), T(i_sig), 'or', 'LineWidth', 2)
                
                while (T(j_zerocrossing) < mean(Stim))
                    j_zerocrossing = j_zerocrossing + 1;
                    if (j_zerocrossing > length(T) / 2)
                        break;
                    end
                    
                end
                
                while ((T(j_sig) < mean(Stim)) && ztest(T(j_sig), baseline_sta, std_baselines, alphaville))
                    
                    j_sig = j_sig + 1;
                    if (j_sig > length(T) / 2)
                        break;
                    end
                    
                end
                
                plot(t(1) - .001 + (j_zerocrossing * .001), T(j_zerocrossing), 'ok', 'LineWidth', 2)
                hold on
                plot(t(1) - .001 + (j_sig * .001), T(j_sig), 'or', 'LineWidth', 2)
                
                if ~ p.vstim
                    
                    peak_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
                    peak_Integration_time_sig = (j_sig - i_sig) * .001;
                    
                elseif p.vstim
                    
                    trough_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
                    trough_Integration_time_sig = (j_sig - i_sig) * .001;
                    
                end
                
                i_zerocrossing_rebound = i_zerocrossing;
                j_zerocrossing_rebound_end = i_zerocrossing_rebound;
                
                while (T(i_zerocrossing_rebound) > mean(Stim))
                    i_zerocrossing_rebound = i_zerocrossing_rebound - 1;
                    if (i_zerocrossing_rebound == 1)
                        break;
                    end
                end
                
                plot(t(1) - .001 + (i_zerocrossing_rebound * .001), T(i_zerocrossing_rebound), 'ok', 'LineWidth', 2)
                plot(t(1) - .001 + ((i_zerocrossing + 1) * .001), T(i_zerocrossing + 1), 'ok', 'LineWidth', 2)
                
                trough_amp_rebound = max(T(i_zerocrossing_rebound:j_zerocrossing_rebound_end));
                
                [H_trough_rebound P_trough_rebound] = ztest(trough_amp_rebound, baseline_sta, std_baselines, alphaville);
                
                if H_trough_rebound
                    
                    if ~ p.vstim
                        trough_Integration_time = (j_zerocrossing_rebound_end - i_zerocrossing_rebound) * .001;
                        
                    elseif p.vstim
                        
                        peak_Integration_time = (j_zerocrossing_rebound_end - i_zerocrossing_rebound) * .001;
                        
                    end
                    
                    trough_location_rebound = find(T(1:length(T) / 2) == trough_amp_rebound);
                    trough_location_time_rebound = abs(ttt(1) + trough_location_rebound * .001 - .001);
                    peak_location_time_rebound = - 10;
                    
                else
                    
                    if ~ p.vstim
                        
                        trough_Integration_time = 0;
                        peak_location_time_rebound = - 10;
                        trough_location_time_rebound = - 10;
                        
                    elseif p.vstim
                        
                        peak_Integration_time = 0;
                        
                    end
                    
                end
                
                if H_trough == 1
                    
                    locs = find(T(1:length(T) / 2) > mean(Stim));
                    i_zerocrossing = (locs(find(locs == trough_location)));
                    j_zerocrossing = (locs(find(locs == trough_location)));
                    i_sig = i_zerocrossing;
                    j_sig = j_zerocrossing;
                    
                    while ((T(i_sig) > mean(Stim)) && ztest(T(i_sig), baseline_sta, std_baselines, alphaville))
                        
                        i_sig = i_sig - 1;
                        if (i_sig == 0)
                            
                            i_sig = 1;
                            break;
                        end
                        
                    end
                    
                    while ((T(j_sig) > mean(Stim)) && ztest(T(j_sig), baseline_sta, std_baselines, alphaville))
                        
                        j_sig = j_sig + 1;
                        if (j_sig > length(T) / 2)
                            break;
                        end
                        
                    end
                    
                    trough_Integration_time_sig = (j_sig - i_sig) * .001;
                    plot(t(1) - .001 + (j_sig * .001), T(j_sig), 'or', 'LineWidth', 2)
                    hold on
                    plot(t(1) - .001 + (i_sig * .001), T(i_sig), 'or', 'LineWidth', 2)
                    
                else
                    
                    trough_Integration_time_sig = 0;
                    
                end
                
            elseif (i_peak < i_trough)
                
                locs = find(T(1:length(T) / 2) > mean(Stim));
                i_zerocrossing = (locs(find(locs == trough_location)));
                j_zerocrossing = (locs(find(locs == trough_location)));
                i_sig = i_zerocrossing;
                j_sig = j_zerocrossing;
                
                while (T(i_zerocrossing) > mean(Stim))
                    
                    i_zerocrossing = i_zerocrossing - 1;
                    
                    if (i_zerocrossing == 0)
                        i_zerocrossing = 1;
                        break;
                        
                    end
                    
                end
                
                while ((T(i_sig) > mean(Stim)) && ztest(T(i_sig), baseline_sta, std_baselines, alphaville))
                    
                    i_sig = i_sig - 1;
                    if (i_sig == 0)
                        i_sig = 1;
                        break;
                    end
                    
                end
                
                plot(t(1) - .001 + (i_zerocrossing * .001), T(i_zerocrossing), 'ok', 'LineWidth', 2)
                hold on
                plot(t(1) - .001 + (i_sig * .001), T(i_sig), 'or', 'LineWidth', 2)
                
                while (T(j_zerocrossing) > mean(Stim))
                    j_zerocrossing = j_zerocrossing + 1;
                    if (j_zerocrossing > length(T) / 2)
                        break;
                    end
                end
                
                while ((T(j_sig) > mean(Stim)) && ztest(T(j_sig), baseline_sta, std_baselines, alphaville))
                    
                    j_sig = j_sig + 1;
                    if (j_sig > length(T) / 2)
                        break;
                    end
                end
                
                plot(t(1) - .001 + (j_zerocrossing * .001), T(j_zerocrossing), 'ok', 'LineWidth', 2)
                hold on
                plot(t(1) - .001 + (j_sig * .001), T(j_sig), 'or', 'LineWidth', 2)
                
                if ~ p.vstim
                    
                    trough_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
                    trough_Integration_time_sig = (j_sig - i_sig) * .001;
                    
                elseif p.vstim
                    
                    peak_Integration_time = (j_zerocrossing - i_zerocrossing) * .001;
                    peak_Integration_time_sig = (j_sig - i_sig) * .001;
                    
                end
                
                i_zerocrossing_rebound = i_zerocrossing;
                j_zerocrossing_rebound_end = i_zerocrossing_rebound;
                
                while (T(i_zerocrossing_rebound) < mean(Stim))
                    i_zerocrossing_rebound = i_zerocrossing_rebound - 1;
                    if (i_zerocrossing_rebound == 1)
                        break;
                    end
                end
                
                plot(t(1) - .001 + (i_zerocrossing_rebound * .001), T(i_zerocrossing_rebound), 'ok', 'LineWidth', 2)
                plot(t(1) - .001 + ((i_zerocrossing + 1) * .001), T(i_zerocrossing + 1), 'ok', 'LineWidth', 2)
                
                peak_amp_rebound = min(T(i_zerocrossing_rebound:j_zerocrossing_rebound_end));
                
                [H_peak_rebound P_peak_rebound] = ztest(peak_amp_rebound, baseline_sta, std_baselines, alphaville);
                
                if H_peak_rebound
                    
                    if ~ p.vstim
                        peak_Integration_time = (j_zerocrossing_rebound_end - i_zerocrossing_rebound) * .001;
                        
                    elseif p.vstim
                        
                        trough_Integration_time = (j_zerocrossing_rebound_end - i_zerocrossing_rebound) * .001;
                        
                    end
                    
                    peak_location_rebound = find(T(1:length(T) / 2) == peak_amp_rebound);
                    peak_location_time_rebound = abs(ttt(1) + peak_location_rebound * .001 - .001);
                    trough_location_time_rebound = - 10;
                    
                else
                    
                    if ~ p.vstim
                        
                        peak_Integration_time = 0;
                        peak_location_time_rebound = - 10;
                        trough_location_time_rebound = - 10;
                        
                    elseif p.vstim
                        
                        trough_Integration_time = 0;
                        
                    end
                    
                end
                
                if H_peak == 1
                    
                    locs = find(T(1:length(T) / 2) < mean(Stim));
                    i_zerocrossing = (locs(find(locs == peak_location)));
                    j_zerocrossing = (locs(find(locs == peak_location)));
                    i_sig = i_zerocrossing;
                    j_sig = j_zerocrossing;
                    
                    while ((T(i_sig) < mean(Stim)) && ztest(T(i_sig), baseline_sta, std_baselines, alphaville))
                        
                        i_sig = i_sig - 1;
                        if (i_sig == 0)
                            
                            i_sig = 1;
                            break;
                        end
                        
                    end
                    
                    while ((T(j_sig) < mean(Stim)) && ztest(T(j_sig), baseline_sta, std_baselines, alphaville))
                        
                        j_sig = j_sig + 1;
                        if (j_sig > length(T) / 2)
                            break;
                        end
                        
                    end
                    
                    peak_Integration_time_sig = (j_sig - i_sig) * .001;
                    plot(t(1) - .001 + (j_sig * .001), T(j_sig), 'or', 'LineWidth', 2)
                    hold on
                    plot(t(1) - .001 + (i_sig * .001), T(i_sig), 'or', 'LineWidth', 2)
                    
                else
                    
                    peak_Integration_time_sig = 0;
                    
                end
                
            elseif (i_peak == i_trough)
                
                peak_Integration_time = 0;
                trough_Integration_time = 0;
                
            end
            
        end
        
    end
    
end

%%

if ~ p.vstim
    
    text(- .8, - 350, 'Peak Location is ')
    text(- .8, - 385, num2str(peak_location_time_rebound))
    text(- .8, - 415, num2str(peak_location_time))
    
    text(- .8, - 450, 'Peak Integration Time is ')
    text(- .8, - 485, num2str(peak_Integration_time))
    text(- .8, - 515, num2str(peak_Integration_time_sig))
    
    text(- .8, - 550, 'Trough Location is ')
    text(- .8, - 585, num2str(trough_location_time_rebound))
    text(- .8, - 615, num2str(trough_location_time))
    
    text(- .8, - 650, 'Trough Integration Time is ')
    text(- .8, - 685, num2str(trough_Integration_time))
    text(- .8, - 715, num2str(trough_Integration_time_sig))
    
elseif p.vstim
    
    text(.4, .9, 'Peak Location is ')
    text(.4, .8, num2str(peak_location_time))
    
    text(.4, .7, 'Peak Integration Time is ')
    text(.4, .6, num2str(peak_Integration_time))
    
    text(.4, .5, 'Trough Location is ')
    text(.4, .4, num2str(trough_location_time))
    
    text(.4, .3, 'Trough Integration Time is ')
    text(.4, .2, num2str(trough_Integration_time))
    
end

if (trough_Integration_time > 0 && peak_Integration_time > 0)
    if peak_location_time < trough_location_time
        
        cell_type = 'ON';
        
    elseif peak_location_time > trough_location_time
        
        cell_type = 'OFF';
        
    end
    
elseif (trough_Integration_time == 0)
    
    cell_type = 'ON';
    
elseif (peak_Integration_time == 0)
    cell_type = 'OFF';
    
end

title(cell_type)

peak_location_time
peak_Integration_time_zero_crossing = peak_Integration_time;
trough_location_time
trough_Integration_time_zero_crossing = trough_Integration_time;

saveas(gcf, [dump, p.cell_id, num2str(p.leave_out), ' ', num2str(p.ii), 'to', num2str(p.jj), '-', 'cubic splined', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), '.fig'], 'fig');
set(gcf, 'PaperPosition', [0 0 20 10]); %x_width=10cm y_width=15cm
saveas(gcf, [dump, p.cell_id, num2str(p.leave_out), ' ', num2str(p.ii), 'to', num2str(p.jj), '-', 'cubic splined', '-BTA is ', num2str(p.cardinal_STA_Only_Burst), '.jpeg'], 'jpeg');

%% p.weighted_burst, p.singleton_spikes
date1 = p.year;
date1(end) = [];

if ~ p.cardinal_STA_Only_Burst
    
    dump = fullfile(p.work_dir, 'Population_Analysis_Alternate\');
    if ~ p.vstim
        file_loc = strcat(dump, p.cell_id, date1, 'single_pulse_analysis = ', num2str(p.single_pulse_activation_correction));
    elseif p.vstim
        file_loc = strcat(dump, p.cell_id, date1);
    end
    
else
    
    dump = fullfile(p.work_dir, 'Population_Analysis_Alternate_BSTA\');
    
    file_loc = strcat(dump, p.cell_id, date1, ' p.weighted_burst = ', num2str(p.weighted_burst), ' p.singleton_spikes = ', num2str(p.singleton_spikes), 'single_pulse_analysis = ', num2str(p.single_pulse_activation_correction));
    
end

if p.vstim
    
    temp_store = peak_T;
    peak_T = trough_T;
    trough_T = temp_store;
    
end

if ~ exist(dump, 'dir'), mkdir(dump); end
save(file_loc, 'peak_location_time', 'trough_location_time', 'cell_type', 'Frequency', 'peak_T', 'trough_T', 'baseline_avg_right', 'baseline_std_right', 'peak_Integration_time_zero_crossing', 'trough_Integration_time_zero_crossing', 'peak_Integration_time_sig', 'trough_Integration_time_sig', 'peak_location_time_rebound', 'trough_location_time_rebound')

end