function pac_segment(subj,block)



if ~strcmp(block(1:3),'Day') && ~strcmp(block(1:3),'All') && ~iscell(block)
    block = {block};
    
else
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



bl = 0.5; %baseline

%if strcmp(task,'RL_baseline')
    ps = 3;
%else
%    ps = 1;   %post-stim
%end


%load all blocks
for iBlock = 1:length(block)
    get_subj_globals(subj,block{iBlock})
    load([dtdir subj '_' block{iBlock} '_theta_gamma_PAC.mat'])
    PAC{iBlock} = ecogPAC;
end

good_elecs = setdiff(banks(1):banks(end),badChannels);

for elec = good_elecs
    
    PAC_cond{elec}.power = [];
    PAC_cond{elec}.bad_trials = [];
    PAC_cond{elec}.cond = [];
    PAC_cond{elec}.value = [];
    PAC_cond{elec}.trl_ind = [];
    
    for iBlock = 1:length(block)
        
        all_block_trials = size(PAC_cond{elec}.power,1);
        get_subj_globals(subj,block{iBlock})
        
        for iTrial = 1:size(allstimtimes,1)
            onset = allstimtimes(iTrial,1);
            PAC_cond{elec}.power(all_block_trials+iTrial,:) = PAC{iBlock}.data(elec,round((onset-bl)*srate):round((onset+ps)*srate));
        end
        
        try
            PAC_cond{elec}.bad_trials = [PAC_cond{elec}.bad_trials PAC{iBlock}.bad_trials{elec}+all_block_trials];
        end
        PAC_cond{elec}.cond = [PAC_cond{elec}.cond ; PAC{iBlock}.cond];
        PAC_cond{elec}.value = [PAC_cond{elec}.value ; PAC{iBlock}.value];
        PAC_cond{elec}.trl_ind = [PAC_cond{elec}.trl_ind ; PAC{iBlock}.trl_ind];
    end
    
end



save([anadir 'PAC_cond_theta_gamma.mat'],'PAC_cond');
end