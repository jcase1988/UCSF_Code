clear
subj = 'EC77';
block = {'B28','B44'};
condition = 10; %0 for stimuli, 1 for value
time_lock = 0; %0 for stim-locked, 1 for consonate-locked
ERP_fit = 0;   %0 for no fitting, 1 for fitting to average ERP
kruns = 10000;

anova_or_regression = 3; %1 for anova, 2 for regression
fig_labels = {'anova','regression','bootstrap'};
fig_label = fig_labels{anova_or_regression};

get_subj_globals(subj,block{1})
good_elecs = banks(1):banks(end);

stimuli = {'aagaa', 'iyfiy','uwshuw'};
valu = {'negative','neutral','positive'};

if time_lock == 1
    %consonant offsets for the three stimuli
    cons_offsets = [round(0.410*srate) round(0.250*srate) round(0.242*srate)];
    lock = 'consonant';
elseif time_lock == 0
    cons_offsets = [0 0 0];
    lock = 'stimulus';
end

if ERP_fit == 0
    fit = 'no_fit';
elseif ERP_fit == 1
    fit = 'fitted';
end


%good_elecs = [243:244 227 211 209 195:196 178:179 149];

bl = 0.5;
ps = 1;

%beg_win = -0.3;  %wrt to time-locked event in secs
%en_win = 0;      %wrt to time-locked event in secs
%beg_win = 0.1;
%en_win = 0.4;
%beg_win = -0.2;
%en_win = -0.15;

beg_win = -0.4;
en_win = 1;
bin_size = 0.1;
dt = 0.02;

%load data from all blocks
for iBlock = 1:length(block)
    get_subj_globals(subj,block{iBlock});
    load([anadir 'HG_cond_power_z-scored.mat'])
    block_data{iBlock} = HG_cond_power;
end

iBin_cnt = 1;
for iBin = beg_win:dt:(en_win-bin_size)
    iBin
        
    %Gather data
    %If multiple blocks
    if iscell(block)
        HG_cat = cell(500,1);
        
        for elec = good_elecs
            
            if ismember(elec,badChannels)
                continue;
            end
            
            for iBlock = 1:length(block)
                get_subj_globals(subj,block{iBlock})
                HG_cond_power = block_data{iBlock};
                                
                if isempty(HG_cond_power{elec})
                    continue;
                end
                
                %turn indices into binary array
                bad_trial_ind = zeros(size(HG_cond_power{elec}.power,1),1);
                bad_trial_ind(HG_cond_power{elec}.bad_trials) = 1;
                HG_cond_power{elec}.bad_trials = bad_trial_ind;
                
                %cancatenate block data (before fitting)
                if isempty(HG_cat{elec})
                    HG_cat{elec} = HG_cond_power{elec};
                else
                    fields = fieldnames(HG_cat{elec})';
                    fields(2,:) = cellfun(@(f) [HG_cat{elec}.(f) ; HG_cond_power{elec}.(f)], fields, 'unif', false);
                    HG_cat{elec} = struct(fields{:});
                end
                
            end
            
            
            if ERP_fit
                [m se trl_m] = fit_ERP(subj,block,HG_cat{elec},iBin,iBin+bin_size,condition,time_lock);
                
                HG_cat{elec}.pow_m = trl_m;
                
            else
                
                if time_lock == 1 %consonant locked
                    
                    for iTrial = 1:size(HG_cat{elec}.power,1)
                        if HG_cat{elec}.cond(iTrial) == 10
                            trl_cond_off = 0;
                        else
                            trl_cond_off = cons_offsets(HG_cat{elec}.cond(iTrial));
                        end
                        tm_win = round(((bl+iBin)*srate)+trl_cond_off):round(((bl+iBin+bin_size)*srate)+trl_cond_off);
                        HG_cat{elec}.pow_m = mean(HG_cat{elec}.power(:,tm_win),2);
                
                    end
                elseif time_lock == 0 %stimulus locked
                    tm_win = round((bl+iBin)*srate:(bl+iBin+bin_size)*srate);
                    HG_cat{elec}.pow_m = mean(HG_cat{elec}.power(:,tm_win),2);
                end
                
                
            end
        end
        
        %If only 1 block
    else
        load([anadir 'HG_cond_power_z-scored.mat'])
        
        %     %turn value indicies into binary array
        %     for elec = good_elecs
        %         value_inds = zeros(size(HG_cond_power{elec}.power,1),1);
        %         value_inds(HG_cond_power{elec}.cond==1) = values(1);
        %         value_inds(HG_cond_power{elec}.cond==2) = values(2);
        %         value_inds(HG_cond_power{elec}.cond==3) = values(3);
        %         value_inds(HG_cond_power{elec}.cond==10) = -999;
        %         HG_cond_power{elec}.values = value_inds;
        %     end
        
        if ERP_fit
            [m se fitted_dat] = fit_ERP(subj,block{iBlock},HG_cond_power{elec},iBin,iBin+bin_size, condition,time_lock);
            HG_cond_power{elec}.power = fitted_dat;
        else
            HG_cat = HG_cond_power;
        end
    end
    


    for elec = good_elecs
        elec
        dat = HG_cat{elec};

        if ismember(elec,badChannels) || size(dat.power,1) < 52
            p_table{iBin_cnt}(1:4,elec) = [elec ; NaN ; NaN ; NaN];
            continue
        end
        
        %remove clicks and first trials
        remove_ind = [find(dat.cond==10) ; find(dat.trl_ind==1)];
        dat.power(remove_ind,:) = [];
        dat.pow_m(remove_ind) = [];
        dat.cond(remove_ind) = [];
        dat.trl_ind(remove_ind) = [];
        dat.bad_trials(remove_ind) = [];
        dat.value(remove_ind) = [];


        clear con_on pow g_stim g_val
        
        pow = dat.pow_m;
        st = dat.cond;
        va = dat.value;

        % Do ANOVA
        if anova_or_regression == 1 
            [p,table,stats,terms] = anovan(pow,{g_stim,g_val},'display','off');
            p_table{iBin_cnt}(1,elec) = elec;
            p_table{iBin_cnt}(2,elec) = p(1);
            p_table{iBin_cnt}(3,elec) = p(2);

            [p,table,stats,terms] = anovan(pow,{g_stim,g_val},'display','off','model','interaction');
            p_table{iBin_cnt}(4,elec) = p(3);
        
        elseif anova_or_regression == 2

            clear st_or va_or X_or
            % order stimulus labels by mean power
            stim_means = [mean(pow(st==1)) mean(pow(st==2)) mean(pow(st==3))];
            [dum i] = sort(stim_means);
            st_or(st==1) = find(i==1);
            st_or(st==2) = find(i==2);
            st_or(st==3) = find(i==3);

            % order value labels by mean power
            value_means = [mean(pow(va==-1)) mean(pow(va==0)) mean(pow(va==1))];
            [dum i] = sort(value_means);
            va_or(va==-1) = find(i==1);
            va_or(va==0) = find(i==2);
            va_or(va==1) = find(i==3);

            % convert to -1,0,1 and creat interaction group
            st_or = st_or-2;
            va_or = va_or-2;
            X_or = st_or .* va_or;

            % n x p matrix for glmfit
            labels = [st_or' va_or' X_or'];
            
            [b,dev,stats] = glmfit(labels,pow);
            
            %add to p_table
            p_table{iBin_cnt}(1,elec) = elec;
            p_table{iBin_cnt}(2:4,elec) = stats.p(2:4);
            
        elseif anova_or_regression == 3
            
            clear samp
            for iStim = 1:3
            
                ind = find(st==iStim);
                st_values = unique(va(ind));
                v1_n = length(find(va(ind)==st_values(1)));
                v2_n = length(find(va(ind)==st_values(2)));
                
                clear samp
                for i = 1:kruns
                    r = randsample(pow(ind),length(ind),true);

                    samp(i) = mean(r(1:v1_n)) -  mean(r(v1_n+1 : length(ind)));
                end
                
                obs1 = pow(st==iStim & va == st_values(1));
                obs2 = pow(st==iStim & va == st_values(2));
                
                obs = mean(obs2) - mean(obs1);
                
                p_table{iBin_cnt}(1,elec) = elec;
                p_table{iBin_cnt}(iStim+1,elec) = length(find(samp<obs))/length(samp);
                
            end                            
            
        end

    end

    for elec = 1:size(p_table{iBin_cnt},2)
        if isnan(p_table{iBin_cnt}(2,elec))
            p_table{iBin_cnt}(2:4,elec) = 1;
        end
    end

    %good_ind = cell2mat(p_table(2:3,:));
    %good_ind = find(good_ind(2,:) < 0.05);
    %p_table(:,good_ind)
    iBin_cnt = iBin_cnt + 1;
end

if ~exist([subjdir [block{:}] '/analysis/'])
    mkdir([subjdir [block{:}] '/analysis/'])
end

save([subjdir [block{:}] '/analysis/' fig_label '_p_table_' [block{:}] '_' lock '_' fit],'p_table');
p_table_heatmap(p_table,subj,block,time_lock,beg_win,en_win,dt,bin_size,fig_label)