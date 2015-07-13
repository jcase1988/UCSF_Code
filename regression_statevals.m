clear
subj = 'EC82';
block = 'B48';
get_subj_globals(subj,block)
load([anadir 'behavior_data.mat'])

stimuli = {'aagaa','iyfiy','uwshuw'};
lock_label = {'per_offered','across_offered','click'};
half_label = {'' '_first_half_only' '_second_half_only'};
locking = 1; %1 for regression per stimID, both offered; 2 for regressions across all stimIDs, both offered; 3 for click for choice
first_second_half_only = 2; %1 first half only, 2 second half only

ps = [0 1]; %window from ps(1) to ps(2)
bl = 0.5;


fig_path = [figdir 'statevals_regression/' lock_label{locking} '_' num2str(ps(1)*1000) '_' num2str(ps(2)*1000) 'ms' half_label{first_second_half_only+1} '/'];
if ~exist(fig_path)
    mkdir(fig_path)
end

load([dtdir subj '_' block '_CAR.mat']);


%statevals represent values AFTER trials...so add [0 0 0] at beginning
statevals_cum = [0 0 0 ; statevals_cum];

%% Setup onsets and Calculate p/beta



if locking == 1 % regression per stimID for both offered
    
    for iStim = 1:3

        titl = stimuli{iStim};
        
        ileft = find(offered(:,1)==iStim);
        iright = find(offered(:,2)==iStim);

        onsets = sort([sample1times(ileft,1) ; sample2times(iright,1)]);
        trl_inds = sort([ileft ; iright]);
        
        %use the 500ms before onset of any sound as the baseline
        for iTrial = 1:length(trl_inds)
            bl_onsets(iTrial) = sample1times(trl_inds(iTrial),1);
        end        

        stateval = statevals_cum(trl_inds,iStim);
        
        if first_second_half_only == 1
            [beta, p] = regression_statevals_helper(subj,ecogCAR,badChannels,bl_onsets(1:round(length(bl_onsets)/2)),onsets(1:round(length(onsets)/2)),bl,ps,per_chan_bad_epochs,stateval, [fig_path stimuli{iStim}],titl);
        elseif first_second_half_only == 2
            [beta, p] = regression_statevals_helper(subj,ecogCAR,badChannels,bl_onsets(round(length(bl_onsets)/2):end),onsets(round(length(onsets)/2):end),bl,ps,per_chan_bad_epochs,stateval, [fig_path stimuli{iStim}],titl);
        else
            [beta, p] = regression_statevals_helper(subj,ecogCAR,badChannels,bl_onsets,onsets,bl,ps,per_chan_bad_epochs,stateval, [fig_path stimuli{iStim}],titl);
        end
        %Plot 
        
        
    end
    
elseif locking == 2 % regression across all stimID for both offered

    titl = 'All-Stim';
    
    onsets = sort([sample1times(:,1) ; sample2times(:,1)]);
    bl_onsets = sort([sample1times(:,1) ; sample1times(:,1)]);
    trl_inds = sort([1:length(sample1times) 1:length(sample1times)]);
    
    cnt = 1;
    for iTrial = 1:2:length(onsets)
        stateval(iTrial) = statevals_cum(cnt,offered(cnt,1));
        stateval(iTrial+1) = statevals_cum(cnt,offered(cnt,2));
        cnt = cnt + 1;
    end
    
    if first_second_half_only == 1        
        [beta, p] = regression_statevals_helper(subj,ecogCAR,badChannels,bl_onsets(1:round(length(bl_onsets)/2)),onsets(1:round(length(onsets)/2)),bl,ps,per_chan_bad_epochs,stateval, [fig_path 'All_stim'],titl);
    elseif first_second_half_only == 2
        [beta, p] = regression_statevals_helper(subj,ecogCAR,badChannels,bl_onsets(round(length(bl_onsets)/2):end),onsets(round(length(onsets)/2):end),bl,ps,per_chan_bad_epochs,stateval, [fig_path 'All_stim'],titl);
    else
       [beta, p] = regression_statevals_helper(subj,ecogCAR,badChannels,bl_onsets,onsets,bl,ps,per_chan_bad_epochs,stateval, [fig_path 'All_stim'],titl);
    end


    
elseif locking == 3 % regression across all stimID for click
        
      
        titl = 'Click';
        
        onsets = clicktimes;
        trl_inds = 1:length(onsets);
        bl_onsets = sample1times(:,1);
        
        for iTrial = 1:length(onsets)
            stateval(iTrial) = statevals_cum(iTrial,offered(iTrial,choice(iTrial))); %find stateval before each click
        end
        
        load([dtdir subj '_' block '_CAR.mat']);
        
        if first_second_half_only == 1
            [beta, p] = regression_statevals_helper(subj,ecogCAR,badChannels,bl_onsets(1:round(length(bl_onsets)/2)),onsets(1:round(length(onsets)/2)),bl,ps,per_chan_bad_epochs,stateval, [fig_path 'Click'],titl);
        elseif first_second_half_only == 2
            [beta, p] = regression_statevals_helper(subj,ecogCAR,badChannels,bl_onsets(round(length(bl_onsets)/2):end),onsets(round(length(onsets)/2):end),bl,ps,per_chan_bad_epochs,stateval, [fig_path 'Click'],titl);
        else
            [beta, p] = regression_statevals_helper(subj,ecogCAR,badChannels,bl_onsets,onsets,bl,ps,per_chan_bad_epochs,stateval, [fig_path 'Click'],titl);
        end
end
    
    
 

