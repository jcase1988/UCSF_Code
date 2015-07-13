function [output_data, output_se, trl_m, output_col] = fit_ERP(subj,block,dat,st_ms,en_ms,condition,time_lock)

%condition = 0 or 1
%0: to plot the three stimuli
%1: to plot the three values that change over blocks

%time_lock = 0 or 1
%0: time-lock to stimulus
%1: time-lock to consonant

%st_ms, en_ms
% in secs, where to start and end with respect to stimulus or consonat (as
% signified by "time_lock"



if iscell(block)
    get_subj_globals(subj,block{1})
else
    get_subj_globals(subj,block);
end


%values = setdiff(unique(dat.value),-999)';
if time_lock == 1
    cons_offsets = [round(0.445*srate) round(0.319*srate) round(0.269*srate)];
elseif time_lock == 0
    cons_offsets = [0 0 0];
end

bl = 0.5 * srate;
labels = {'aagaa','iyfiy','uwshuw'};

st_sa = round(st_ms*srate);
en_sa = round(en_ms*srate);

% first realign data with respect to the three stimuli,
% regardless of single or multiple blocks
[output_data_stim, trl_m] = fit_ERP_helper(dat,bl,cons_offsets,st_sa,en_sa);

%if rows of output should correspond to AAGAA-B1, AAGAA-B2, IYFIY-B1, IYFIY-B2, UWSHUW-B1, UWSHUW-B2
if condition == 3
    value_cols = {'-b','-r','-g'};
    
    cond_block_inds = [1 1 ; 1 2 ; 2 1 ; 2 2 ; 3 1 ; 3 2];
    
    for iComp = 1:size(cond_block_inds,1)

        dat_val = output_data_stim(dat.cond==cond_block_inds(iComp,1) & dat.iBlock==cond_block_inds(iComp,2),:);
        output_data(iComp,:) = nanmean(dat_val,1);
        output_se(iComp,:) = nanstd(dat_val,[],1)/sqrt(size(dat_val,1));
        
        if cond_block_inds(iComp,2) == 1
            output_col{iComp} = '-k';
        else
            output_col{iComp} = value_cols{unique(dat.value(dat.cond==cond_block_inds(iComp,1) & dat.iBlock==2))+2};
        end
    end

   
    
    
%if rows of output should correspond to AAGAA-NEG, AAGAA-NEUTRAL, IYFIY-NEG, IYFIY-POS, UWSHUW - NEUTRAL, UWSHUW - POS
elseif condition == 2
    
    %codes for aagaa-neg, aagaa-neutral, etc.
    cond_val_inds = unique([dat.cond dat.value],'rows');
    cond_val_inds(cond_val_inds(:,1) == 10,:) = [];
    value_cols = {'-b','-r','-g'};
    
    for iComp = 1:size(cond_val_inds,1);    
        dat_val = output_data_stim(dat.cond==cond_val_inds(iComp,1) & dat.value==cond_val_inds(iComp,2) & dat.bad_trials==0 & dat.iBlock~=1,:); %iBlock = 1 for pre-learning day1
        output_data(iComp,:) = nanmean(dat_val,1);
        output_se(iComp,:) = nanstd(dat_val,[],1)/sqrt(size(dat_val,1));
        output_col{iComp} = value_cols{cond_val_inds(iComp,2)+2};
        
    end

    %last three rows will be aagaa,iyfiy,uwshuw for pre-learning day 1    
    for iCond = 1:3
        dat_val = output_data_stim(dat.cond==iCond & dat.bad_trials==0 & dat.iBlock==1,:); %iBlock = 1 for pre-learning day1
        output_data(iComp+iCond,:) = nanmean(dat_val,1);
        output_se(iComp+iCond,:) = nanstd(dat_val,[],1)/sqrt(size(dat_val,1));
    end
    
%if rows of output should correspond to NEGATIVE, NEUTRAL, POSITIVE
elseif condition == 1
    if length(unique(values)) > 1 %if values are not all neutral, make sure output is sorted from negative-netural-positive
        [sorted_values] = sort(values);
        
        for iVal = 1:length(sorted_values)
            
            dat_val = output_data_stim(dat.value==sorted_values(iVal) & dat.bad_trials==0,:); %will not index clicks
            
            output_data(iVal,:) = nanmean(dat_val,1); %average across blocks
            output_se(iVal,:) = nanstd(dat_val,1)/sqrt(size(dat_val,1));
        end
    elseif unique(values) == 0
        dat_val = output_data_stim(dat.value==0 & dat.bad_trials==0,:);
        output_data(2,:) = nanmean(dat_val,1);
        output_se(2,:) = nanstd(dat_val,1)/sqrt(size(dat_val,1));
    end
    
%if rows of output should correspond to AAGAA, IYFIY, UWSHUW
elseif condition == 0
    for Stim = 1:3
        
        dat_stim = output_data_stim(dat.cond==Stim & dat.bad_trials==0,:); %will not index clicks
        
        output_data(Stim,:) = nanmean(dat_stim,1); %average across blocks
        output_se(Stim,:) = nanstd(dat_stim,1)/sqrt(size(dat_stim,1));
    end
else
    output_data = [];
    output_se = [];
end
end


function [output_data, trl_m] = fit_ERP_helper(dat,bl,cons_offsets,st_sa,en_sa)

for iStim = 1:3
    tm_win = 1+(bl+cons_offsets(iStim)+st_sa):(bl+cons_offsets(iStim)+en_sa)+1;
    
    
    inds = find(dat.cond == iStim)';
    
    
    val_pow = dat.power(inds,:);
    fitted_pow = val_pow(:,tm_win); %just to start off
    
    lags = zeros(size(val_pow,1),1);
    for iteration = 1:5
        avg_ERP = mean(fitted_pow,1);
        
        for i = 1:size(fitted_pow,1)
            [c l] = xcorr(fitted_pow(i,:),avg_ERP,40);
            [a max_corr_ind] = max(c);
            
            if (lags(i) + l(max_corr_ind)+tm_win(end)) <= size(val_pow,2) && (lags(i) + l(max_corr_ind) + tm_win(1)) >= 1            
                lags(i) = lags(i) + l(max_corr_ind);
            end
        end
        
        for iTrial = 1:size(val_pow,1)
            fitted_pow(iTrial,:) = val_pow(iTrial,tm_win+lags(iTrial));
        end

    end
    
    cnt = 1;
    for i = inds
        trl_m(i,1) = mean(fitted_pow(cnt,:));
        output_data(i,:) = fitted_pow(cnt,:);
        cnt = cnt + 1;
    end
    
    
end
% bad_trials = find(trl_m==0);
% trl_m(bad_trials) = [];
% output_data(bad_trials,:) = [];
end

