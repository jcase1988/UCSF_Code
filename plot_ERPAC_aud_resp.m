function plot_ERPAC_aud_resp(subj,block,condition,elecs)

%condition = 0 or 1 or 2
%0: to plot the three stimuli
%1: to plot the three values that change over blocks
%-1: to plot no condition

if ~iscell(block)
    get_subj_globals(subj,block)
    block = {block};
else
    get_subj_globals(subj,block{1})
    figdir = [subjdir [block{:}] '/figures/'];
end

fig_path = [figdir 'line_graph_theta_gamma_ERPAC/'];

if condition == 0
    fig_path = [fig_path 'stim/'];
elseif condition == 1
    fig_path = [fig_path 'value/'];
    
elseif condition == -1
    fig_path = [fig_path 'no_conds/'];
end

if ~exist(fig_path)
    mkdir(fig_path)
end

% phase_freq = [2.5 5];
% amp_freq = [34 130];
phase_freq = [4 7];
amp_freq = [70 150];


srate = 400;
bl = 0.5;
ps = 2;
timeWindow = [-bl ps]*1000;

start_time_window = -bl*1000;
plot_jump = 250;
tm_st=1;
jm = 0.25*srate;
tm_en = (bl+ps)*srate;

for z = 1:length([tm_st:jm:tm_en])
    plot_str{z} = round(start_time_window+(z-1)*plot_jump);
end

data = [];
eventsA = [];
eventsB = [];
eventsC = [];
block_len = 0;

for iBlock = 1:length(block)
    get_subj_globals(subj,block{iBlock})
    load([dtdir subj '_' block{iBlock} '_CAR.mat'])
    
    if condition == 0
        eventsA = [eventsA round(allstimtimes(find(stimID==1),1)' * srate + block_len)];
        eventsB = [eventsB round(allstimtimes(find(stimID==2),1)' * srate + block_len)];
        eventsC = [eventsC round(allstimtimes(find(stimID==3),1)' * srate + block_len)];        
        
    elseif condition == 1
        eventsA = [eventsA round(allstimtimes(find(stimID==find(values==-1)),1)' * srate + block_len)];
        eventsB = [eventsB round(allstimtimes(find(stimID==find(values==0)),1)' * srate + block_len)];
        eventsC = [eventsC round(allstimtimes(find(stimID==find(values==1)),1)' * srate + block_len)];
    elseif condition == -1
        events = allstimtimes(:,1)'*srate + block_len;
    end
        
    %concatenate data across blocks
    data = [data ecogCAR.data];
    block_len = block_len + size(ecogCAR.data,2);
end
    
if ~exist('elecs')
    elecs = 1:size(data,1) ;
end

for elec = elecs
    
    if condition ~= -1
        [pacA, pacB, pacC] = erpac_corr3(data(elec,:), srate, eventsA, eventsB, eventsC, timeWindow, phase_freq(1), phase_freq(2), amp_freq(1), amp_freq(2));

        %figure;
        figure('visible','off');
        h = plot(1:length(pacA),band_pass(pacA,srate,0,7),'-b'); set(h,'LineWidth',2)
        hold on;
        h = plot(1:length(pacB),band_pass(pacB,srate,0,7),'-r'); set(h,'LineWidth',2)
        hold on;
        h = plot(1:length(pacC),band_pass(pacC,srate,0,7),'-g'); set(h,'LineWidth',2)
        
        cond_label = 'conds';
    else
        pac = erpac_corr(data(elec,:), srate, events, timeWindow, phase_freq(1), phase_freq(2), amp_freq(1), amp_freq(2));

        %figure;
        figure('visible','off');
        h = plot(1:length(pac),band_pass(pac,srate,0,20),'-b'); set(h,'LineWidth',2)
        cond_label = 'no_conds';
    end
        
    set(gca, 'XTick', [tm_st:jm:tm_en], 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');
    
    if condition == 0
        legend('aagaa','iyfiy','uwshuw')
    elseif condition == 1
        legend('negative','neutral','positive')
    end
    xlabel('ms');
    ylabel('PAC');
    xlim([tm_st tm_en])
    title(['e' num2str(elec)])
    
    saveas(gcf, [fig_path 'ERPAC_' cond_label '_e' num2str(elec) '.jpg'], 'jpg')
    
    close;
         
end
    

