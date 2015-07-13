function pac_within_elec(subj,block,coupling)

%coupling: 1 for theta-HG coupling

if coupling == 1
    phase_freq = [4 7];
    power_freq = [70 150];
    
    phase_lab = 'theta';
    power_lab = 'gamma';
end

get_subj_globals(subj,block)
load([dtdir subj '_' block '_CAR.mat'])
ecogCAR = notf(ecogCAR);
ecogPAC = rmfield(ecogCAR,'data');

bl = 0.5; %baseline

if strcmp(task,'RL_baseline')
    ps = 3;
else
    ps = 1;   %post-stim
end


%if LEARNING block
if exist('sample1times1times')
    onsets = [sample1times(:,1) ; sample2times(:,1) ; feedbacktimes(:,1)]*srate;
    offsets = [sample1times(:,2) ; sample2times(:,2) ; feedbacktimes(:,2)]*srate;
%if BASELINE block
else
    onsets = round(allstimtimes(:,1)*srate);
    offsets = round(allstimtimes(:,2)*srate);
end

ecogPAC.cond = [];
ecogPAC.value = [];
ecogPAC.trl_ind = [];

good_elecs = setdiff(banks(1):banks(end),badChannels);

for iTrial = 1:length(onsets)
    ecogPAC.cond = [ecogPAC.cond ; stimID(iTrial)];
    
    %if one of the three stimuli
    if ismember(stimID(iTrial),[1 2 3])
        ecogPAC.value(iTrial,1) = values(stimID(iTrial));
    %if click
    else
        ecogPAC.value(iTrial,1) = -999;
    end   
    
    ecogPAC.trl_ind(iTrial,1) = iTrial;
    
end



for elec = good_elecs
    ecogPAC.bad_trials{elec} = [];
    for iTrial = 1:length(onsets)
        
        bl_onset = floor(onsets(iTrial)-(bl*srate));
        win = floor(bl_onset:onsets(iTrial)+(ps*srate));   
        
        %identify bad trials
        for iEpoch = 1:size(per_chan_bad_epochs{elec},1)
            beg = per_chan_bad_epochs{elec}(iEpoch,1)*srate;
            en = per_chan_bad_epochs{elec}(iEpoch,2)*srate;

            if (beg > win(1) && beg < win(end)) || (en > win(1) && en < win(end))
                ecogPAC.bad_trials{elec} = [ecogPAC.bad_trials{elec} iTrial];
            end
        end
    end
end


for iChan = 1:size(ecogCAR.data,1)
    iChan
    signal = ecogCAR.data(iChan,:);
    
    % Extract low frequency analytic phase
    phase_freq_phase = eegfilt(signal, srate, phase_freq(1), []);
    phase_freq_phase = eegfilt(phase_freq_phase, srate, [], phase_freq(2));
    phase_freq_phase = angle(hilbert(phase_freq_phase));
    
    % Extract low frequency analytic phase of high frequency analytic amplitude
    power_freq_phase = eegfilt(signal, srate, power_freq(1), []);
    power_freq_phase = eegfilt(power_freq_phase, srate, [], power_freq(2));
    power_freq_phase = abs(hilbert(power_freq_phase));
    power_freq_phase = eegfilt(power_freq_phase, srate, phase_freq(1), []);
    power_freq_phase = eegfilt(power_freq_phase, srate, [], phase_freq(2));
    power_freq_phase = angle(hilbert(power_freq_phase));
    
    % Calculate PAC
    ecogPAC.data(iChan,:) = exp(1i * (phase_freq_phase - power_freq_phase));
end




save([dtdir subj '_' block '_' phase_lab '_' power_lab '_PAC.mat'],'ecogPAC')

end