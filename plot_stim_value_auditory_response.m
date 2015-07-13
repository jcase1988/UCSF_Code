function plot_stim_value_auditory_response(subj,block,time_lock,ERP_fit,elecs_to_plot)

%block = {Pre-learning-Day1, etc.}

%time_lock = 0 or 1
%0: time-lock to stimulus
%1: time-lock to consonant 

%ERP_fit = 0 or 1
%0: do not do average response fitting
%1: fit to average response within each class


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


fig_path = [figdir 'stim_value_interaction/'];

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

dat_label = 'Z-score';

%consonant offsets for the three stimuli
cons_offsets = [round(0.445*srate) round(0.319*srate) round(0.269*srate)];

bl = 0.5; %baseline

if strcmp(task,'RL_baseline')
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
    tm_en = round( (ps-0.5)*srate);
    plot_jump = 250;
    jm = round(plot_jump./1000*srate);

    tm_en = tm_en + 40;
end


%XTick labels
for z = 1:length([tm_st:jm:tm_en])
    plot_str{z} = round(start_time_window+(z-1)*plot_jump);
end

%Gather data

HG = cell(500,1);

%HG only power


for iBlock = 1:length(block)
    get_subj_globals(subj,block{iBlock})
    
    load([anadir 'HG_cond_power_z-scored.mat'])    
    HG_cond_power_cell{iBlock} = HG_cond_power;
end

if exist('elecs_to_plot')
    good_elecs = elecs_to_plot;
else
    good_elecs = setdiff(banks(1):banks(end),badChannels);
end

for elec = good_elecs
    
    for iBlock = 1:length(block)
        
        if ~isempty(HG_cond_power_cell{iBlock}{elec})
                %turn indices into binary array
                bad_trial_ind = zeros(size(HG_cond_power_cell{iBlock}{elec}.power,1),1);
                bad_trial_ind(HG_cond_power_cell{iBlock}{elec}.bad_trials) = 1;
                HG_cond_power_cell{iBlock}{elec}.bad_trials = bad_trial_ind;


            if iBlock == 1 || isempty(HG{elec})
                HG{elec} = HG_cond_power_cell{iBlock}{elec};
                HG{elec}.iBlock = iBlock*ones(size(HG_cond_power_cell{iBlock}{elec}.power,1),1);
            else
                HG_cond_power_cell{iBlock}{elec}.iBlock = ones(size(HG_cond_power_cell{iBlock}{elec}.power,1),1)*iBlock;
                fields = fieldnames(HG{elec})';
                fields(2,:) = cellfun(@(f) [HG{elec}.(f) ; HG_cond_power_cell{iBlock}{elec}.(f)], fields, 'unif', false);
                HG{elec} = struct(fields{:});
            end
        end
        
    end
end



for elec = good_elecs
    dat = HG{elec};
    
    if isempty(dat) %|| length(unique(HG{elec}.iBlock))~=3
        continue;
    end
    
    if ERP_fit == 0 %do not fit to average ERP

        %Reject bad trials
        %A = aagaa, B = iyfiy, C = uwshuw
        if condition == 0
            good_trials_A = setdiff(find(dat.cond==1 & dat.trl_ind~=1 & dat.iBlock ~=1),find(dat.bad_trials));
            good_trials_B = setdiff(find(dat.cond==2 & dat.trl_ind~=1),find(dat.bad_trials));
            good_trials_C = setdiff(find(dat.cond==3 & dat.trl_ind~=1),find(dat.bad_trials));

        %A = negative, B = neutral, C = positive values
        elseif condition == 1
            good_trials_A = setdiff(find(dat.value==-1 & dat.trl_ind~=1 & dat.iBlock ~=1),find(dat.bad_trials));
            good_trials_B = setdiff(find(dat.value==0 & dat.trl_ind~=1 & dat.iBlock ~=1),find(dat.bad_trials));
            good_trials_C = setdiff(find(dat.value==1 & dat.trl_ind~=1 & dat.iBlock ~=1),find(dat.bad_trials));
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
        
        
        
        figure;
        
        %AAGAA for NEG and NEUTRAL
        subplot(3,1,1)
        h = boundedline(tm_st:tm_en,band_pass(power_mean(1,:),srate,0,7),band_pass(power_se(1,:),srate,0,7),'alpha','transparency',0.15,'-b'); set(h,'LineWidth',2);
        hold on;
        h = boundedline(tm_st:tm_en,band_pass(power_mean(2,:),srate,0,7),band_pass(power_se(2,:),srate,0,7),'alpha','transparency',0.15,'-r'); set(h,'LineWidth',2);
        hold on;
        h = boundedline(tm_st:tm_en,band_pass(power_mean(7,:),srate,0,7),band_pass(power_se(7,:),srate,0,7),'alpha','transparency',0.15,'-k'); set(h,'LineWidth',2);
        title([num2str(elec) ' - AAGAA'])
     %   legend('Negative','Neutral')
        xlabel('ms');
        ylabel('Z-Score'); 
        xlim([tm_st tm_en]) 
        set(gca, 'XTick', [tm_st:jm:tm_en], 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');    
        
        %IYFIY for NEG and POS
        subplot(3,1,2)
        h = boundedline(tm_st:tm_en,band_pass(power_mean(3,:),srate,0,7),band_pass(power_se(3,:),srate,0,7),'alpha','transparency',0.15,'-b'); set(h,'LineWidth',2);
        hold on;
        h = boundedline(tm_st:tm_en,band_pass(power_mean(4,:),srate,0,7),band_pass(power_se(4,:),srate,0,7),'alpha','transparency',0.15,'-g'); set(h,'LineWidth',2);
        hold on;
        h = boundedline(tm_st:tm_en,band_pass(power_mean(2,:),srate,0,7),band_pass(power_se(8,:),srate,0,7),'alpha','transparency',0.15,'-k'); set(h,'LineWidth',2);
        title([num2str(elec) ' - IYFIY'])
      %  legend('Negative','Positive')
        xlabel('ms');
        ylabel('Z-Score'); 
        xlim([tm_st tm_en]) 
        set(gca, 'XTick', [tm_st:jm:tm_en], 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');    
        
        %UWHSUW for NEUTRAL and POS
        subplot(3,1,3)
        h = boundedline(tm_st:tm_en,band_pass(power_mean(5,:),srate,0,7),band_pass(power_se(5,:),srate,0,7),'alpha','transparency',0.15,'-r'); set(h,'LineWidth',2);
        hold on;
        h = boundedline(tm_st:tm_en,band_pass(power_mean(6,:),srate,0,7),band_pass(power_se(6,:),srate,0,7),'alpha','transparency',0.15,'-g'); set(h,'LineWidth',2);
        hold on;
        h = boundedline(tm_st:tm_en,band_pass(power_mean(2,:),srate,0,7),band_pass(power_se(9,:),srate,0,7),'alpha','transparency',0.15,'-k'); set(h,'LineWidth',2);
        title([num2str(elec) ' - UWSHUW'])
    %    legend('Neutral','Positive')   
        xlabel('ms'); 
        ylabel('Z-Score'); 
        xlim([tm_st tm_en]) 
        set(gca, 'XTick', [tm_st:jm:tm_en], 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');    
        
    
        
        saveas(gcf, [fig_path dat_label '_conds_e' num2str(elec) '.png'], 'png')
        close;
        
        
        
        
        
     elseif ERP_fit == 1
%         
        [power_mean, power_se, m, col] = fit_ERP(subj,block,dat,-0.5,ps,2,time_lock);

    % Plot Z-scored signal for each condition    
    
  
  
    
        figure;
        
        %AAGAA for NEG and NEUTRAL
        subplot(3,1,1)
        h = boundedline(tm_st:tm_en,band_pass(power_mean(1,:),srate,0,7),band_pass(power_se(1,:),srate,0,7),'alpha','transparency',0.15,col{1}); set(h,'LineWidth',2);
        hold on;
        h = boundedline(tm_st:tm_en,band_pass(power_mean(2,:),srate,0,7),band_pass(power_se(2,:),srate,0,7),'alpha','transparency',0.15,col{2}); set(h,'LineWidth',2);
        hold on;
        h = boundedline(tm_st:tm_en,band_pass(power_mean(7,:),srate,0,7),band_pass(power_se(7,:),srate,0,7),'alpha','transparency',0.15,'-k'); set(h,'LineWidth',2);
        title([num2str(elec) ' - AAGAA'])
     %   legend('Negative','Neutral')
        xlabel('ms');
        ylabel('Z-Score'); 
        xlim([tm_st tm_en]) 
        set(gca, 'XTick', [tm_st:jm:tm_en], 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');    
        
        %IYFIY for NEG and POS
        subplot(3,1,2)
        h = boundedline(tm_st:tm_en,band_pass(power_mean(3,:),srate,0,7),band_pass(power_se(3,:),srate,0,7),'alpha','transparency',0.15,col{3}); set(h,'LineWidth',2);
        hold on;
        h = boundedline(tm_st:tm_en,band_pass(power_mean(4,:),srate,0,7),band_pass(power_se(4,:),srate,0,7),'alpha','transparency',0.15,col{4}); set(h,'LineWidth',2);
        hold on;
        h = boundedline(tm_st:tm_en,band_pass(power_mean(8,:),srate,0,7),band_pass(power_se(8,:),srate,0,7),'alpha','transparency',0.15,'-k'); set(h,'LineWidth',2);
        title([num2str(elec) ' - IYFIY'])
      %  legend('Negative','Positive')
        xlabel('ms');
        ylabel('Z-Score'); 
        xlim([tm_st tm_en]) 
        set(gca, 'XTick', [tm_st:jm:tm_en], 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');    
        
        %UWHSUW for NEUTRAL and POS
        subplot(3,1,3)
        h = boundedline(tm_st:tm_en,band_pass(power_mean(5,:),srate,0,7),band_pass(power_se(5,:),srate,0,7),'alpha','transparency',0.15,col{5}); set(h,'LineWidth',2);
        hold on;
        h = boundedline(tm_st:tm_en,band_pass(power_mean(6,:),srate,0,7),band_pass(power_se(6,:),srate,0,7),'alpha','transparency',0.15,col{6}); set(h,'LineWidth',2);
        hold on;
        h = boundedline(tm_st:tm_en,band_pass(power_mean(9,:),srate,0,7),band_pass(power_se(9,:),srate,0,7),'alpha','transparency',0.15,'-k'); set(h,'LineWidth',2);
        title([num2str(elec) ' - UWSHUW'])
    %    legend('Neutral','Positive')   
        xlabel('ms'); 
        ylabel('Z-Score'); 
        xlim([tm_st tm_en]) 
        set(gca, 'XTick', [tm_st:jm:tm_en], 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');    
        
    
        
        saveas(gcf, [fig_path dat_label '_conds_e' num2str(elec) '.png'], 'png')
        close;
    end
end
