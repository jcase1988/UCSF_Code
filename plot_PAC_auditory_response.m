function plot_PAC_auditory_response(subj,block,condition,time_lock,coupling)

%condition = 0 or 1 or 2
%0: to plot the three stimuli
%1: to plot the three values that change over blocks
%2: to plot one subplot for each stimuli with lines for each instance of value 

%time_lock = 0 or 1
%0: time-lock to stimulus
%1: time-lock to consonant 

%Coupling = 0 or 1
%0: PLV HG power plot 
%1: PLV theta-HG phase-amplitude coupling


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

if coupling == 0
    fig_path = [figdir 'line_graph_aud_resp/'];
    dat_label = 'Z-score';
elseif coupling == 1
    fig_path = [figdir 'line_graph_theta_gamma_PAC_binned/'];
    dat_label = 'PAC';
end

if condition == 0
    fig_path = [fig_path 'stim/'];
elseif condition == 1
    fig_path = [fig_path 'value/'];
end

if ~exist(fig_path)
    mkdir(fig_path)
end



if iscell(block)
    for iBlock = 1:length(block)
        get_subj_globals(subj,block{iBlock})
        load([anadir 'PAC_cond_theta_gamma.mat'])
        PAC{iBlock} = PAC_cond;
    end
    
    for elec = 1:length(PAC{1})
        PAC_cat{elec}.power = [];
        PAC_cat{elec}.cond = [];
        PAC_cat{elec}.value = [];
        PAC_cat{elec}.trl_ind = [];
        PAC_cat{elec}.bad_trials = [];
        
        trl_cnt = 0;
        for iBlock = 1:length(block)
            if ~isempty(PAC{iBlock}{elec})
                PAC_cat{elec}.power = [PAC_cat{elec}.power ; PAC{iBlock}{elec}.power];
                PAC_cat{elec}.cond = [PAC_cat{elec}.cond ; PAC{iBlock}{elec}.cond];
                PAC_cat{elec}.value = [PAC_cat{elec}.value ; PAC{iBlock}{elec}.value];
                PAC_cat{elec}.trl_ind = [PAC_cat{elec}.trl_ind ; PAC{iBlock}{elec}.trl_ind];
                
                PAC_cat{elec}.bad_trials = [PAC_cat{elec}.bad_trials  PAC{iBlock}{elec}.bad_trials+trl_cnt];
                trl_cnt = trl_cnt + size(PAC_cat{elec}.power,1);
            end
        end
        if isempty(PAC_cat{elec}.power)
            PAC_cat(elec) = [];
        end
    end
        
    PAC_cond = PAC_cat;
        
else
    load([anadir 'PAC_cond_theta_gamma.mat'])
end
% 
% if time_lock == 0
% 
%     start_time_window = -bl*1000;
%     tm_st = round( start_time_window ./1000 *srate);
%     tm_en = round( ps*srate);
%     plot_jump = 250;
%     jm = round(plot_jump./1000*srate);
% 
% elseif time_lock == 1
%     
%     start_time_window = -bl*1000;
%     tm_st = round( start_time_window ./1000 *srate);
%     tm_en = round( (ps-0.5)*srate);
%     plot_jump = 250;
%     jm = round(plot_jump./1000*srate);
% 
%     tm_en = tm_en + 40;
% end



win_size = 0.05 * srate; %50 ms
win_dt = 0.01 * srate; %10 ms
bins = 1:win_dt:(size(PAC_cond{elec}.power,2)-win_size);

xticks_labels = -500:250:2500;
%XTick labels
for z = 1:length(xticks_labels)
    [a i] = min(abs((xticks_labels(z)/1000) - ((bins/srate)-0.5)));
    xticks(z) = i;
end



for elec = 1:length(PAC_cond)
    dat = PAC_cond{elec};
    
    if isempty(dat)
        continue;
    end
    
    

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
    
   
    
    
    for iBin = 1:length(bins)
        samps = dat.power(good_trials_A,bins(iBin):bins(iBin)+win_size);        
        pac_A(iBin) = mean(abs(mean(samps,1)));
        angles_A{iBin} = angle(samps);
        
        samps = dat.power(good_trials_B,bins(iBin):bins(iBin)+win_size);        
        pac_B(iBin) = mean(abs(mean(samps,1)));        
        angles_B{iBin} = angle(samps);
        
        samps = dat.power(good_trials_C,bins(iBin):bins(iBin)+win_size);        
        pac_C(iBin) = mean(abs(mean(samps,1)));
        angles_C{iBin} = angle(samps);     
    

    end

    % Plot PAC time course   
    
    h = figure;
    h = plot(1:length(pac_A),pac_A,'-b'); set(h,'LineWidth',2)
    hold on;
    h = plot(1:length(pac_B),pac_B,'-r'); set(h,'LineWidth',2)
    hold on;
    h = plot(1:length(pac_C),pac_C,'-g'); set(h,'LineWidth',2)
    set(gca, 'XTick', xticks, 'XTickLabel', xticks_labels, 'XTickMode', 'manual', 'Layer', 'top');
    
    if condition == 0
        legend('aagaa','iyfiy','uwshuw')
    elseif condition == 1
        legend('negative','neutral','positive')
    end
    xlabel('ms');
    ylabel('PAC');
    xlim([1 (2*srate)/win_dt])
    title(['e' num2str(elec)])        
        
    saveas(gcf, [fig_path 'PAC_conds_e' num2str(elec) '.png'], 'png')
    close;
    
    
    
    % Plot PAC rose and line plot
    
%     
%     for iBin = 1:length(bins)
%         h = figure;
%         plot(1:length(pac_A),pac_A,'-b'); hold on; plot([iBin iBin],[0 max(pac_A)],'-b'); set(h,'LineWidth',2)
%         set(gca, 'XTick', xticks, 'XTickLabel', xticks_labels, 'XTickMode', 'manual', 'Layer', 'top');
%         title('Negative')
%         
% %         scrollsubplot(2,3,(iBin-1)*6+2)    
%         hold on;
%         plot(1:length(pac_B),pac_B,'-r'); hold on; plot([iBin iBin],[0 max(pac_B)],'-b'); set(h,'LineWidth',2)
%         set(gca, 'XTick', xticks, 'XTickLabel', xticks_labels, 'XTickMode', 'manual', 'Layer', 'top');
%         title('Neutral')
%         
% %         scrollsubplot(2,3,(iBin-1)*6+3)        
%         plot(1:length(pac_C),pac_C,'-g'); hold on; plot([iBin iBin],[0 max(pac_C)],'-b'); set(h,'LineWidth',2)
%         set(gca, 'XTick', xticks, 'XTickLabel', xticks_labels, 'XTickMode', 'manual', 'Layer', 'top');
%         title('Positive')
%             
% %         xlabel('ms');
% %         ylabel('PAC');
% %         xlim([1 (2*srate)/win_dt])
% %         title(['e' num2str(elec)])     
%         
% %         scrollsubplot(2,3,(iBin-1)*6+4)
% %         rose(angles_A{iBin}(:))
% %         scrollsubplot(2,3,(iBin-1)*6+5)
% %         rose(angles_B{iBin}(:))
% %         scrollsubplot(2,3,(iBin-1)*6+6)
% %         rose(angles_C{iBin}(:))
%     end
        
        
%     saveas(gcf, [fig_path 'PAC_conds_e' num2str(elec) '.png'], 'png')
%     close;
    
end
