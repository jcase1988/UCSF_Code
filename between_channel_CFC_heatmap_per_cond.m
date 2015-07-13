function between_channel_CFC_heatmap_per_cond(subj,block,condition,coupling,e_low,e_high)

% 1. Computes instanteous PLV then computes the correlation (beta) between
% value [negative (-1), neutral (0), positive (1)] and PLV (averaged across
% trials) for a particular time bin. Each electrode pair will have 1 beta.

% 2. Plot a heat map Betas in an adjacency matrix
%
%

%condition: 0 for stimulus
%           1 for value




%coupling: 1 for theta-gamma coupling
%          2 for theta-theta coupling
%          3 for alpha-alpha coupling
%          4 for low gamma - low gamma


%e_low/e_high: electrodes for phase-providing frequency
%                    and  for power-providing frequency


if ~iscell(block)
    get_subj_globals(subj,block)
    
    % Reassign paths if "All" or "Day" Variable
    if strcmp(block(1:3),'All')
        anadir = anadir_all;
        figdir = figdir_all;
    elseif strcmp(block(1:3),'Day')
        anadir = anadir_day;
        figdir = figdir_day;
    end
    %if "block" is a list of blocks to concatenate
else
    get_subj_globals(subj,block{1});
    figdir = [subjdir [block{:}] '/figures/'];
end

if coupling == 1
    fig_path = [figdir 'conditions_btw_channel_PAC_theta_gamma/'];
    dat_label = 'PAC_theta_gamma';
    freq1 = [4 7];
    freq2 = [70 150];
    title_label = 'Theta Gamma PAC';
elseif coupling == 2
    fig_path = [figdir 'conditions_btw_channel_PLV_theta_theta/'];
    dat_label = 'PLV_theta_theta';
    freq = [4 7];
    title_label = 'Theta Theta PLV';
elseif coupling == 3
    fig_path = [figdir 'conditions_btw_channel_PLV_alpha_alpha/'];
    dat_label = 'PLV_alpha_alpha';
    freq = [8 12];
    title_label = 'Alpha Alpha PLV';
elseif coupling == 4
    fig_path = [figdir 'conditions_btw_channel_PLV_lowgamma_lowgamma/'];
    dat_label = 'PLV_lowgamma_lowgamma';
    freq = [30 70];
    title_label = 'Low Gamma Low Gamma PLV';
end

%setup labels for the beginning of each strip/depth
if strcmp(subj,'EC77')
    bank_labels = {};
    bank_labels(1:256) = {'Grid'};
    bank_labels(257:260) = {'AOF'};
    bank_labels(261:264) = {'POF'};
    bank_labels(265:270) = {'ITG'};
    bank_labels(271:276) = {'Pos-Temp'};
    bank_labels(277:280) = {'AST'};
    bank_labels(281:284) = {'MST'};
    bank_labels(285:288) = {'PST'};
    bank_labels(289:292) = {'Insula'};
    bank_labels(293:296) = {'sACC'};
    bank_labels(297:300) = {'iACC'};
    bank_labels(301:304) = {'Amyg'};
    bank_labels(305:308) = {'Hippo'};
end


if coupling ~= 1
    e_high = e_low;
    
    prev_plotted = {};
    e_label_x = {};
    e_label_y = {};
    for e = e_high
        if ismember(bank_labels(e),prev_plotted)
            e_label_x = [e_label_x ' '];
            e_label_y = [e_label_y {num2str(e)}];
        else
            e_label_x = [e_label_x {bank_labels{e}(1:2)}];
            e_label_y = [e_label_y {[bank_labels{e} ' ' num2str(e)]}];
            
            prev_plotted = [prev_plotted {bank_labels{e}}];
        end
    end
end

if condition == 0 
    fig_path = [fig_path 'stim/'];
    condition_labels = {'aagaa','iyfiy','uwshuw'};
elseif condition == 1
    fig_path = [fig_path 'value/'];
    condition_labels = {'negative','neutral','positive'};
end

%fig_path = [fig_path 'e' num2str(e_low(1)) '_e' num2str(e_low(end)) '/'];

if ~exist(fig_path)
    mkdir(fig_path)
end




bl = 0.5;
ps = 2;

if ~iscell(block)
    block = {block};
end



for e1 = e_low
    for e2 = e_high
        PLV{e1,e2}.PLV = [];
        PLV{e1,e2}.cond = [];
        PLV{e1,e2}.value = [];
        PLV{e1,e2}.trl_ind = [];
        PLV{e1,e2}.bad_trials = [];
    end
end

all_trial_cnt = 0;
for iBlock = 1:length(block)
    clear value dat_low dat_high
    get_subj_globals(subj,block{iBlock})
    load([dtdir subj '_' block{iBlock} '_CAR.mat'])
    
    %Phase-amplitude coupling
    if coupling == 1
        
        %calculate phase locking between electrodes
        dat_low(e_low,:) = eegfilt(ecogCAR.data(e_low,:),srate,freq1(1),[]); %high-pass
        dat_low(e_low,:) = eegfilt(dat_low(e_low,:),srate,[],freq1(2)); %low-pass
        dat_low(e_low,:) = angle(hilbert(dat_low(e_low,:))); %instanteous phase
        
        dat_high(e_high,:) = eegfilt(ecogCAR.data(e_high,:),srate,freq2(1),[]);
        dat_high(e_high,:) = eegfilt(dat_high(e_high,:),srate,[],freq2(2));
        dat_high(e_high,:) = abs(hilbert(dat_high(e_high,:)));
        dat_high(e_high,:) = eegfilt(dat_high(e_high,:),srate,freq1(1),[]);
        dat_high(e_high,:) = eegfilt(dat_high(e_high,:),srate,[],freq1(2));
        dat_high(e_high,:) = angle(hilbert(dat_high(e_high,:)));
        
        %PLV
    elseif coupling == 2 || coupling == 3 || coupling == 4
        
        %calculate phase locking between electrodes
        
        for e1 = e_low
            dat_low(e1,:) = eegfilt(ecogCAR.data(e1,:),srate,freq(1),[]); %high-pass
            dat_low(e1,:) = eegfilt(dat_low(e1,:),srate,[],freq(2)); %low-pass
            dat_low(e1,:) = angle(hilbert(dat_low(e1,:))); %instanteous phase
        end
        
        dat_high = dat_low;
        
    end
    
    %calculate value (1,2,3,-999) for each trial
    for i = 1:length(stimID)
        try
            value(i) = values(stimID(i))+2;
        catch
            value(i) = -999;
        end
    end
    
    for e1 = e_low
        if coupling ~= 1
            e_high = e_low(e_low>e1);
            if isempty(e_high)
                e_high = e_low;
                continue
            end
        end
        for e2 = e_high
            %calculate PLV
            dat = exp(1i * (dat_low(e1,:) - dat_high(e2,:)));
            
            PLV{e1,e2}.trl_ind = [PLV{e1,e2}.trl_ind ; (1:length(stimID))'];
            
            for iTrial = 1:size(allstimtimes,1)
                wind = round((allstimtimes(iTrial,1)-bl)*srate):round((allstimtimes(iTrial,1)+ps)*srate);
                
                PLV{e1,e2}.PLV = [PLV{e1,e2}.PLV ; dat(wind)];
                PLV{e1,e2}.cond = [PLV{e1,e2}.cond ; stimID(iTrial)];
                PLV{e1,e2}.value = [PLV{e1,e2}.value ; value(iTrial)];
                
                %identify bad trials
                bad_epochs = [per_chan_bad_epochs{e1} ; per_chan_bad_epochs{e2}];
                for iEpoch = 1:size(bad_epochs,1)
                    beg = bad_epochs(iEpoch,1)*srate;
                    en = bad_epochs(iEpoch,2)*srate;
                    
                    if (beg > wind(1) && beg < wind(end)) || (en > wind(1) && en < wind(end))
                        PLV{e1,e2}.bad_trials = [PLV{e1,e2}.bad_trials iTrial+all_trial_cnt];
                    end
                end
                
                
            end
        end
    end
    
    all_trial_cnt = all_trial_cnt + length(stimID);
end


%plot PLV

win_size = 0.25 * srate; %250 ms
win_dt = 0.25 * srate; %250 ms
bins = 1:win_dt:(size(PLV{e_low(1),e_high(2)}.PLV,2)-win_size);

bin_times = -500:250:1750;

for iCond = 1:length(condition_labels)

    for iBin = 1:length(bins)
        clear y_PLV
        e1_cnt = 1;

        for e1 = e_low
            if coupling ~= 1
                e_high = e_low(e_low>e1);
                e2_cnt = e1_cnt+1;

                if isempty(e_high)
                    e_high = e_low;
                    continue
                end
            end
            for e2 = e_high
                dat = PLV{e1,e2};

                if condition == 0
                    good_trials = setdiff(find(dat.cond==iCond),dat.bad_trials);

                elseif condition == 1
                    good_trials = setdiff(find(dat.value==iCond),dat.bad_trials);

                end

                samps = mean(dat.PLV(good_trials,bins(iBin):bins(iBin)+win_size),1);
                y_PLV(e1_cnt,e2_cnt) = mean(abs(samps));
                y_PLV(e2_cnt,e1_cnt) = y_PLV(e1_cnt,e2_cnt);
                e2_cnt = e2_cnt + 1;
            end
            e1_cnt = e1_cnt + 1;
        end
        
        
        figure('visible','off');
        heatmap(y_PLV,e_label_x,e_label_y); colorbar;
        caxis([0 0.5])
        title([num2str(bin_times(iBin)) ' ms']);
        saveas(gcf, [fig_path condition_labels{iCond} '_e' num2str(e_low(1)) '-e' num2str(e_low(end)) '_' num2str(bin_times(iBin)) 'ms.jpg'], 'jpg')
        close;
        
        
    end


end




