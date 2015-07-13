function [ecogSeg] = segment_trial(ecogCAR,bl,ps,allstimtimes)


%bl = 0.5; %baseline in seconds
%ps = 3; % trial offset in seconds
srate = ecogCAR.sampFreq;
events = allstimtimes;

nchan = size(ecogCAR.data,1);
nsamp = round((bl + ps) * srate)+1;
ntrial = size(events,1);

ecogSeg = zeros(nchan,nsamp,ntrial);

for iTrial = 1:size(allstimtimes,1)
    onset = events(iTrial,1);
    offset = events(iTrial,2);
    
    trl_win = round(((-bl + onset) * srate)):round(((onset + ps) * srate));
    if length(trl_win) > nsamp
        trl_win = trl_win(1:nsamp);
    end
    ecogSeg(:,:,iTrial) = ecogCAR.data(:,trl_win);
end
end
