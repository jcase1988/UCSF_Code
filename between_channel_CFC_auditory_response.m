function between_channel_CFC_auditory_response(subj,block,condition,coupling,e_low,e_high)

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
    fig_path = [figdir 'line_graph_btw_channel_PAC_theta_gamma/'];
    dat_label = 'PAC_theta_gamma';
    freq1 = [4 7];
    freq2 = [70 150];
    title_label = 'Theta Gamma PAC';
elseif coupling == 2
    fig_path = [figdir 'line_graph_btw_channel_PLV_theta_theta/'];
    dat_label = 'PLV_theta_theta';
    freq = [4 7];
    title_label = 'Theta Theta PLV';
elseif coupling == 3
    fig_path = [figdir 'line_graph_btw_channel_PLV_alpha_alpha/'];
    dat_label = 'PLV_alpha_alpha';
    freq = [8 12];
    title_label = 'Alpha Alpha PLV';
elseif coupling == 4
    fig_path = [figdir 'line_graph_btw_channel_PLV_lowgamma_lowgamma/'];
    dat_label = 'PLV_lowgamma_lowgamma';
    freq = [30 70];
    title_label = 'Low Gamma Low Gamma PLV';
end

if coupling ~= 1
    e_high = e_low;
end

if condition == 0
    fig_path = [fig_path 'stim/'];
elseif condition == 1
    fig_path = [fig_path 'value/'];
end

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
        
        for e1 = e_low
            %calculate phase locking between electrodes
            dat_low(e1,:) = eegfilt(ecogCAR.data(e1,:),srate,freq1(1),[]); %high-pass
            dat_low(e1,:) = eegfilt(dat_low(e1,:),srate,[],freq1(2)); %low-pass
            dat_low(e1,:) = angle(hilbert(dat_low(e1,:))); %instanteous phase
        end
        
        for e2 = e_high
            dat_high(e2,:) = eegfilt(ecogCAR.data(e2,:),srate,freq2(1),[]);
            dat_high(e2,:) = eegfilt(dat_high(e2,:),srate,[],freq2(2));
            dat_high(e2,:) = abs(hilbert(dat_high(e2,:)));
            dat_high(e2,:) = eegfilt(dat_high(e2,:),srate,freq1(1),[]);
            dat_high(e2,:) = eegfilt(dat_high(e2,:),srate,[],freq1(2));
            dat_high(e2,:) = angle(hilbert(dat_high(e2,:)));
        end
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

win_size = 0.05 * srate; %50 ms
win_dt = 0.01 * srate; %10 ms
bins = 1:win_dt:(size(PLV{e_low(1),e_high(2)}.PLV,2)-win_size);

xticks_labels = -500:250:2000;
%XTick labels
for z = 1:length(xticks_labels)
    [a i] = min(abs((xticks_labels(z)/1000) - ((bins/srate)-0.5)));
    xticks(z) = i;
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
        dat = PLV{e1,e2};
        
        if condition == 0
            good_trials_A = setdiff(find(dat.cond==1),dat.bad_trials);
            good_trials_B = setdiff(find(dat.cond==2),dat.bad_trials);
            good_trials_C = setdiff(find(dat.cond==3),dat.bad_trials);
        elseif condition == 1
            good_trials_A = setdiff(find(dat.value==1),dat.bad_trials);
            good_trials_B = setdiff(find(dat.value==2),dat.bad_trials);
            good_trials_C = setdiff(find(dat.value==3),dat.bad_trials);
        end
        
        for iBin = 1:length(bins)
            samps = mean(dat.PLV(good_trials_A,bins(iBin):bins(iBin)+win_size),1);
            plv_A(iBin) = mean(abs(samps));
            
            samps = mean(dat.PLV(good_trials_B,bins(iBin):bins(iBin)+win_size),1);
            plv_B(iBin) = mean(abs(samps));
            
            samps = mean(dat.PLV(good_trials_C,bins(iBin):bins(iBin)+win_size),1);
            plv_C(iBin) = mean(abs(samps));
        end
        
        figure('visible','off');
        h = plot(1:length(plv_A),plv_A,'-b'); set(h,'LineWidth',2)
        hold on;
        h = plot(1:length(plv_B),plv_B,'-r'); set(h,'LineWidth',2)
        hold on;
        h = plot(1:length(plv_C),plv_C,'-g'); set(h,'LineWidth',2)
        set(gca, 'XTick', xticks, 'XTickLabel', xticks_labels, 'XTickMode', 'manual', 'Layer', 'top');
        
        if condition == 0
            legend('aagaa','iyfiy','uwshuw')
        elseif condition == 1
            legend('negative','neutral','positive')
        end
        xlabel('ms');
        ylabel('PAC');
        xlim([1 length(plv_A)])
        title([title_label ' e' num2str(e1) ' - e' num2str(e2)])
        
        saveas(gcf, [fig_path dat_label '_e' num2str(e1) '-_e' num2str(e2) '.jpg'], 'jpg')
        close;
        
    end
end
end

