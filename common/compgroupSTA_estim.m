function [STA, stimulus_matrix, len, NSP, inner, Stm_matrix] = compgroupSTA_estim(Stm, Spikes, tKerLen, STA, NSP, inner, ceiling)
%function [STA, stimulus_matrix, len, NSP, inner] = compgroupSTA_estim(Stm, Spikes, tKerLen, STA, NSP, inner, ceiling)
%
% COMPutes GROUP (multiple repeats of the stimulus) Spike Triggered
%  Average for Electrical STIMulation data.
%
%  Inputs:  
%    Stm = stimulus (2D matrx: first index is time bin index, second is space)
%    (units are the amplitude units of the stimulus, e.g. luminance or
%    voltage)
%    Spikes = column vector of spike counts in each time bin (length of 'Stm') 
%    tKerLen = length of time window to be used for plotting the STA.
%    Expresed in terms of bins. Hence for a frequency of 50Hz stimulation
%    tKerLen of 50 represents a window of 1s for which the STA is
%    calculated. Since frequency of 50Hz is 20ms, and hence 50 bins is
%    1second
%    STA = Vector of stimuli preceding spikes in given bin multiplied by the number of
%    spikes in that given bin
%    NSP = Number of Spikes across all trials used for the final STA
%    calculation
%    inner = total number of stimuli vectors used for STA calculation
%    across trials
%    ceiling = length of stimulus

% This code is mirrored compgroupSTA.m for visual stimuli.

%CODE

Spinds = find(Spikes(1:end));%Find which positions in the spike histogram had spikes
length(Spinds);
location=find(Spinds<=tKerLen);%Since I want to store tKerLen frames preceding each spike, I need to exclude those spikes which ocured within the 
% first tKerLen frames in a stimulus trial. So I find the location of all the spikes that occured
% before the first tKerLen frames

%% Finds the number of spikes that happened within the last tKerLen from the end of the stimulus.


upperlocation=find(Spinds>(ceiling-tKerLen));
uppercut=length(upperlocation);

if isempty(location)
    
  loc=0;
    
else
    
    loc=length(location);% Find the last spike that occured before the first tKerLen frames and find its position
end

%% Create a stimulus matrix of the appropriate length (by removing the spikes that occured too early or too late i.e which occured within tKerLen frames of the beginning or end of 
%% stimulus block

len=length(Spinds)-loc-uppercut; % number of stimuli vectors to be included in STA calculation in this trial
stimulus_matrix=zeros(len,2*tKerLen); % creation of STA vector. The STA vector is tKerLen before 0s and tKerLen after 0s

inner=inner+len; % number of stimuli vectors to be included in STA calculation across all trials

[slen,swid]= size(Stm);

nsp = 0; % number of spikes in each trial
Stm_matrix = [];


%% If there are no spikes it returns this message and skips the STA calculation step
counter=1;
if isempty(Spinds)
    fprintf(1, 'Not enough spikes to compute STA - compSTA\n');
   
else

    %% For loop that actually calculates the STA
    for i = 1:length(Spinds)


        if((Spinds(i))>(tKerLen)&&(Spinds(i))<(ceiling-tKerLen)) % If the spike occurs after the first tKerLen frames in a trial or before the last tKerLen frames in a trial include the stimuli for the STA calculation

        STA = STA + (Stm(:,(Spinds(i))-tKerLen+1:(Spinds(i))+tKerLen).*Spikes(Spinds(i)))'; % The STA vector adds up all the stimuli occuring before a spike. The stimuli are multiplied by the number of spikes they cause before they
        % are used for the STA calculation

        stimulus_matrix(counter,:)=Stm(:,(Spinds(i))-tKerLen+1:(Spinds(i))+tKerLen)'; % Stores the stimuli used for the STA calculation

        
        counter=counter+1;
        
        nsp=nsp+Spikes(Spinds(i)); % Update the number of spikes used for the STA calculation

        
        stm_matrix = zeros(Spikes(Spinds(i)),tKerLen*2);

        
         for t = 1:Spikes(Spinds(i))
             
             stm_matrix(t,:)=Stm(:,(Spinds(i))-tKerLen+1:(Spinds(i))+tKerLen)';
             
             
         end
         
         Stm_matrix = [Stm_matrix; stm_matrix];

        else
            
        end
         
    end
    
    
    
end
%% Finds total number of spikes across trials which have been used to calculate the STA
NSP=nsp+NSP;