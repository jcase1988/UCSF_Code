function z_score_auditory_response(subj,block)

get_subj_globals(subj,block)

if ~exist(anadir)
    mkdir(anadir)
end

if ~iscell(block)
    blocks = {block};
end


get_subj_globals(subj,blocks{1}) %to get elecs
HG_cond_power = cell(banks(end),1);
trial_total_cnt = 0; %for indexing bad trials
 
for iBlock = 1:length(blocks)
    keep subj block coupling blocks iBlock trial_total_cnt HG_cond_power anadir_all
    get_subj_globals(subj,blocks{iBlock})
    
    
    load([dtdir subj '_' blocks{iBlock} '_CAR.mat'])
    HG = hgf(ecogCAR);
    label = 'power';

    

    bl = 0.5; %baseline
    
    if strcmp(lower(task),'rl_baseline')
        ps = 3;
    else
        ps = 1;   %post-stim
    end
    
    good_elecs = setdiff(banks(1):banks(end),badChannels);

    %if LEARNING block
    if exist('sample1times1times')
        onsets = [sample1times(:,1) ; sample2times(:,1) ; feedbacktimes(:,1)]*srate;
        offsets = [sample1times(:,2) ; sample2times(:,2) ; feedbacktimes(:,2)]*srate;
      %  stimID(stimID~=10) %exclude click onsets

    %if BASELINE block
    else
        onsets = round(allstimtimes(:,1)*srate);
        offsets = round(allstimtimes(:,2)*srate);

        %exclude click onsets
     %   onsets(stimID == 10) = [];
     %   offsets(stimID == 10) = [];
     %   stimID(stimID~=10); %exclude click onsets
    end

    

    for elec = good_elecs
        clear band
        band = HG.data(elec,:);

        if isempty(HG_cond_power{elec})
            HG_cond_power{elec}.power = [];
            HG_cond_power{elec}.cond = [];
            HG_cond_power{elec}.value = [];
            HG_cond_power{elec}.trl_ind = [];
            HG_cond_power{elec}.bad_trials = [];
        end
        
        
        for iTrial = 1:length(onsets)

            %baseline gamma signal
            bl_onset = floor(onsets(iTrial)-(bl*srate));
            bl_win = floor(bl_onset:onsets(iTrial));
            win = floor(bl_onset:onsets(iTrial)+(ps*srate));            

            %z-transform
            bl_m = mean(band(bl_win));
            bl_sd = std(band(bl_win));
            zed = (band(win)-bl_m)./bl_sd;            

            %assign z-scored signal, condition, and trial index
            HG_cond_power{elec}.power = [HG_cond_power{elec}.power ; zed];
            HG_cond_power{elec}.cond = [HG_cond_power{elec}.cond ; stimID(iTrial)];
            HG_cond_power{elec}.trl_ind = [HG_cond_power{elec}.trl_ind ; iTrial];
            
            %if one of the three stimuli
            if ismember(stimID(iTrial),[1 2 3])
                HG_cond_power{elec}.value(iTrial,1) = values(stimID(iTrial));
            %if click
            else
                HG_cond_power{elec}.value(iTrial,1) = -999;
            end
           
            %identify bad trials
            for iEpoch = 1:size(per_chan_bad_epochs{elec},1)
                beg = per_chan_bad_epochs{elec}(iEpoch,1)*srate;
                en = per_chan_bad_epochs{elec}(iEpoch,2)*srate;                
                
                if (beg > win(1) && beg < win(end)) || (en > win(1) && en < win(end))
                    HG_cond_power{elec}.bad_trials = [HG_cond_power{elec}.bad_trials trial_total_cnt+iTrial];
                end                
            end                                    
        end
      
    end
    
    trial_total_cnt = trial_total_cnt + length(onsets);
end

if length(block)>2
    
    if strcmp(block(1:3),'Day')
        anadir = anadir_day;
        
        if ~exist(anadir)
            mkdir(anadir)
        end
    elseif strcmp(block(1:3),'All')
        anadir = anadir_all;
        
        if ~exist(anadir)
            mkdir(anadir)
        end
    end
end

save([anadir 'HG_cond_' label '_z-scored.mat'],'HG_cond_power');
return




