function plot_power_auditory_response(subj,block,condition,time_lock,ERP_fit,coupling,elecs_to_plot)

%condition = 0 or 1 or 2
%0: to plot the three stimuli
%1: to plot the three values that change over blocks
%2: to plot one subplot for each stimuli with lines for each instance of value 

%time_lock = 0 or 1
%0: time-lock to stimulus
%1: time-lock to consonant 

%ERP_fit = 0 or 1
%0: do not do average response fitting
%1: fit to average response within each class

%Coupling = 0 or 1
%0: HG power plot
%1: theta-HG phase-amplitude coupling

if ~iscell(block)
    get_subj_globals(subj,block)

    % Reassign paths if "All" or "Day" Variable
    try
        if strcmp(block(1:3),'All')
            anadir = anadir_all;
            figdir = figdir_all;
        elseif strcmp(block(1:3),'Day')
            anadir = anadir_day;
            figdir = figdir_day;
        end
    end
%if "block" is a list of blocks to concatenate
else 
    get_subj_globals(subj,block{1});
    figdir = [subjdir [block{:}] '/figures/'];
end

if coupling == 0
    fig_path = [figdir 'line_graph_aud_resp/'];
    dat_label = 'Z-score';
elseif coupling == 1
    fig_path = [figdir 'line_graph_theta_gamma_PAC/'];
    dat_label = 'PAC';
end


if condition == 0
    fig_path = [fig_path 'stim/'];
elseif condition == 1
    fig_path = [fig_path 'value/'];
elseif condition == 2
    fig_path = [fig_path 'stim_value_interaction/'];
end

if time_lock == 0
    fig_path = [fig_path 'stim-locked'];
elseif time_lock == 1
    fig_path = [fig_path 'consonant'];
end

if ERP_fit == 1
    fig_path = [fig_path '_fitted/'];
elseif ERP_fit == 0
    fig_path = [fig_path '/'];
end

if ~exist(fig_path)
    mkdir(fig_path)
end

%consonant offsets for the three stimuli
cons_offsets =  [round(0.410*srate) round(0.250*srate) round(0.242*srate)];


bl = 0.5; %baseline

if strcmp(lower(task),'rl_baseline')
    ps = 1.5;
else
    ps = 1;   %post-stim
end


if time_lock == 0

    start_time_window = -bl*1000;
    tm_st = round( start_time_window ./1000 *srate);
    tm_en = round( ps*srate);
    plot_jump = 250;
    jm = round(plot_jump./1000*srate);

elseif time_lock == 1
    
    start_time_window = -bl*1000;
    tm_st = round( start_time_window ./1000 *srate);
    tm_en = round( (ps)*srate);
    plot_jump = 250;
    jm = round(plot_jump./1000*srate);

end


%XTick labels
for z = 1:length([tm_st:jm:tm_en])
    plot_str{z} = round(start_time_window+(z-1)*plot_jump);
end

%Gather data
if iscell(block)
    HG = cell(500,1);
    
    %HG only power
    if coupling == 0

        for iBlock = 1:length(block)
            get_subj_globals(subj,block{iBlock})
            
            load([anadir 'HG_cond_power_z-scored.mat'])
            
            if exist('elecs_to_plot')
                good_elecs = elecs_to_plot;
            else                
                good_elecs = setdiff(banks(1):banks(end),badChannels);
            end
            
            for elec = good_elecs

                if isempty(HG_cond_power{elec})                    
                    continue;
                end
                
                %turn indices into binary array
                bad_trial_ind = zeros(size(HG_cond_power{elec}.power,1),1);
                bad_trial_ind(HG_cond_power{elec}.bad_trials) = 1;
                HG_cond_power{elec}.bad_trials = bad_trial_ind;

                if isempty(HG{elec})
                    HG{elec} = HG_cond_power{elec};
                else
                    fields = fieldnames(HG{elec})';
                    fields(2,:) = cellfun(@(f) [HG{elec}.(f) ; HG_cond_power{elec}.(f)], fields, 'unif', false);
                    HG{elec} = struct(fields{:});
                end

                %HG{elec} = [HG{elec} HG_cond_power{elec}];

            end
        end
        
    elseif coupling == 1
    
        %load all blocks
        for iBlock = 1:length(block)
            get_subj_globals(subj,block{iBlock})
            load([dtdir subj '_' block{iBlock} '_theta_gamma_PAC.mat'])            
            PAC{iBlock} = ecogPAC;
        end
        
        if exist('elecs_to_plot')
            good_elecs = elecs_to_plot;
        else
            good_elecs = setdiff(banks(1):banks(end),badChannels);
        end
        
        for elec = good_elecs
            
            HG{elec}.power = [];
            HG{elec}.bad_trials = [];
            HG{elec}.cond = [];
            HG{elec}.value = [];
            HG{elec}.trl_ind = [];
            
            for iBlock = 1:length(block)
                
                all_block_trials = size(HG{elec}.power,1);
                get_subj_globals(subj,block{iBlock})
                                
                for iTrial = 1:size(allstimtimes,1)                    
                    onset = allstimtimes(iTrial,1);
                    HG{elec}.power(all_block_trials+iTrial,:) = PAC{iBlock}.data(elec,round((onset-0.5)*srate):round((onset+2)*srate));                    
                end
                
                try
                    HG{elec}.bad_trials = [HG{elec}.bad_trials PAC{iBlock}.bad_trials{elec}+all_block_trials];
                end
                HG{elec}.cond = [HG{elec}.cond ; PAC{iBlock}.cond];
                HG{elec}.value = [HG{elec}.value ; PAC{iBlock}.value];
                HG{elec}.trl_ind = [HG{elec}.trl_ind ; PAC{iBlock}.trl_ind];
            end
%         
%             if iBlock == 1
%                 HG{elec} = HG_cond_power{elec};
%             else
%                 fields = fieldnames(HG{elec})';
%                 fields(2,:) = cellfun(@(f) [HG{elec}.(f) ; HG_cond_power{elec}.(f)], fields, 'unif', false);
%                 HG{elec} = struct(fields{:});
%             end
        end
    end
    
else
    load([anadir 'HG_cond_power_z-scored.mat'])
    good_elecs = setdiff(banks(1):banks(end),badChannels);
    
    

    
    %turn value indicies into binary array
%     for elec = good_elecs
%         value_inds = zeros(size(HG_cond_power{elec}.power,1),1);
%         value_inds(HG_cond_power{elec}.cond==1) = values(1);
%         value_inds(HG_cond_power{elec}.cond==2) = values(2);
%         value_inds(HG_cond_power{elec}.cond==3) = values(3);
%         value_inds(HG_cond_ower{elec}.cond==10) = -999;
%         HG_cond_power{elec}.values = value_inds;
%     end
    
    HG = HG_cond_power;
end


for elec = good_elecs
    dat = HG{elec};
    
    if isempty(dat)
        continue;
    end
    
    if ~iscell(block)
        %turn indices into binary array
        bad_trial_ind = zeros(size(HG_cond_power{elec}.power,1),1);
        bad_trial_ind(HG_cond_power{elec}.bad_trials) = 1;
        dat.bad_trials = bad_trial_ind;
    end
    

    
    if ERP_fit == 0 %do not fit to average ERP

        %Reject bad trials
        %A = aagaa, B = iyfiy, C = uwshuw
        if condition == 0
            good_trials_A = setdiff(find(dat.cond==1 & dat.trl_ind~=1),find(dat.bad_trials));
            good_trials_B = setdiff(find(dat.cond==2 & dat.trl_ind~=1),find(dat.bad_trials));
            good_trials_C = setdiff(find(dat.cond==3 & dat.trl_ind~=1),find(dat.bad_trials));

        %A = negative, B = neutral, C = positive values
        elseif condition == 1
            good_trials_A = setdiff(find(dat.value==-1 & dat.trl_ind~=1),find(dat.bad_trials));
            good_trials_B = setdiff(find(dat.value==0 & dat.trl_ind~=1),find(dat.bad_trials));
            good_trials_C = setdiff(find(dat.value==1 & dat.trl_ind~=1),find(dat.bad_trials));
        end

        %Retrieve HG envelope   
        %if time-locked to stimulus
        if time_lock == 0
            power_A = dat.power(good_trials_A,1:srate*(bl+ps)+1);
            power_B = dat.power(good_trials_B,1:srate*(bl+ps)+1);
            power_C = dat.power(good_trials_C,1:srate*(bl+ps)+1);
        elseif time_lock == 1
            power_A = []; power_B = []; power_C = [];
            for iTrial = good_trials_A'
                con_on = bl*srate + cons_offsets(dat.cond(iTrial));
                power_A = [power_A ; dat.power(iTrial, round((con_on-srate*0.5):(con_on+srate*ps)))];
            end

            for iTrial = good_trials_B'
                con_on = bl*srate + cons_offsets(dat.cond(iTrial));
                power_B = [power_B ; dat.power(iTrial, round((con_on-srate*0.5):(con_on+srate*ps)))];
            end

            for iTrial = good_trials_C'
                con_on = bl*srate + cons_offsets(dat.cond(iTrial));
                power_C = [power_C ; dat.power(iTrial, round((con_on-srate*0.5):(con_on+srate*ps)))];
            end
        end

        if coupling == 0
            %Calculate mean and SE over time
            power_A_se = nanstd(power_A,1)/sqrt(size(power_A,1));
            power_B_se = nanstd(power_B,1)/sqrt(size(power_B,1));
            power_C_se = nanstd(power_C,1)/sqrt(size(power_C,1));
            power_A = nanmean(power_A,1);
            power_B = nanmean(power_B,1);
            power_C = nanmean(power_C,1);  
        else
            power_A_se = 0;
            power_B_se = 0;
            power_C_se = 0;
            power_A = abs(nanmean(power_A,1));
            power_B = abs(nanmean(power_B,1));
            power_C = abs(nanmean(power_C,1)); 
        end

        if isnan(power_A(1))
            continue;
        end
    elseif ERP_fit == 1
        
        [power_mean, power_se, m] = fit_ERP(subj,block,dat,-0.5,ps,condition,time_lock);
        
        if condition == 0 || condition == 1

            power_A = power_mean(1,:); %negative / aagaa
            power_B = power_mean(2,:); %neutral  / iyfiy
            power_C = power_mean(3,:); %positive / uwshuw

            power_A_se = power_se(1,:); %negative / aagaa
            power_B_se = power_se(2,:); %neutral  / iyfiy
            power_C_se = power_se(3,:); %positive / uwshuw
        end
            
       %if condition == 2, do seperate plot for aagaa-neg, aagaa-neutral, etc. 
        
    end

    % Plot Z-scored signal for each condition    
    
    
    figure;
    h = boundedline(tm_st:tm_en,band_pass(power_A,srate,0,7),band_pass(power_A_se,srate,0,7),'alpha','transparency',0.15,'-b'); set(h,'LineWidth',2);
    hold on;
    h = boundedline(tm_st:tm_en,band_pass(power_B,srate,0,7),band_pass(power_B_se,srate,0,7),'alpha','transparency',0.15,'-r'); set(h,'LineWidth',2);
    hold on;
    h = boundedline(tm_st:tm_en,band_pass(power_C,srate,0,7),band_pass(power_C_se,srate,0,7),'alpha','transparency',0.15,'-g'); set(h,'LineWidth',2);
    set(gca, 'XTick', [tm_st:jm:tm_en], 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');
    
    if condition == 0
        legend('aagaa','iyfiy','uwshuw')
    elseif condition == 1
        legend('negative','neutral','positive')
    end
    xlabel('ms');
    ylabel('Z-Score');
    xlim([tm_st tm_en])
    title(['e' num2str(elec)])
    
    
        
    saveas(gcf, [fig_path dat_label '_conds_e' num2str(elec) '.png'], 'png')
    close;
end
