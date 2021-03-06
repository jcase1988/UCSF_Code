function plot_power_grid(subj,block)
get_subj_globals(subj,block)

% Reassign paths if "All" or "Day" Variable
if length(block)>2
    if strcmp(block(1:3),'All')
        anadir = anadir_all;
        figdir = figdir_all;
    elseif strcmp(block(1:3),'Day')
        anadir = anadir_day;
        figdir = figdir_day;
    end
end

fig_path = [figdir 'power_grid_aud_resp/'];

if ~exist(fig_path)
    mkdir(fig_path)
end

good_elecs = setdiff(banks(1):banks(end),badChannels);

if strcmp(task,'RL_baseline')
    bl = 0.5; %baseline
    ps = 2;   %post-stim
elseif strcmp(task,'RL_learning')
    bl = 0.5;
    ps = 1;
end

start_time_window = -bl*1000;
tm_st = round( start_time_window ./1000 *srate);
tm_en = round( ps*srate);
plot_jump = 250;
jm = round(plot_jump./1000*srate);

ScSz=[1 1 2400 1260];
pos_grid=[1 1 2400 1260];


if strcmp(subj,'EC71')
    elecs = {1:64, 65:128};
    elec_locs = {[1:8 ; 9:16 ; 17:24 ; 25:32 ; 33:40 ; 41:48 ; 49:56 ; 57:64], [1:8 ; 9:16 ; 17:24 ; 25:32 ; 33:40 ; 41:48 ; 49:56 ; 57:64]};
    y_spacing = 35;
    dimens = [8 8];
    grid_labels = {'OFC_grid','lateral_grid'};
    
else
    elecs = {1:256};
    elec_locs = {[1:16 ; 17:32 ; 33:48 ; 49:64 ; 65:80 ; 81:96 ; 97:112 ; 113:128 ; 129:144 ; 145:160 ; 161:176 ; 177:192 ; 193:208 ; 209:224 ; 225:240 ; 241:256]};
    y_spacing = 30;
    dimens = [16 16];
    grid_labels = {'grid'};
    
end




for iGrid = 1:length(elecs)
    
    
    elec_locs{iGrid} = flipud(flipud(elec_locs{iGrid})')';
    if strcmp(subj,'EC82') || strcmp(subj,'EC71') %if left grid
        elec_locs{iGrid} = fliplr(elec_locs{iGrid});
    end
    elec_locs{iGrid} = reshape(elec_locs{iGrid},1,length(elecs{iGrid}));
    
    
    gridPos = DivideScreen(dimens(1),dimens(2), ScSz, 20, y_spacing, pos_grid); % Change first two variables to your first grid/strip's dimensions.
    
    fgrid = figure('visible','off');
    set(fgrid, 'Position', ScSz,'color', [1 1 1], 'MenuBar', 'figure');
    
    
    cnt = 1;
    
    load([anadir 'HG_cond_power_z-scored.mat'])
    es = elecs{iGrid};%intersect(elecs,good_elecs);
    for iElec = 1:length(es)
        clear HG
        dat = HG_cond_power{es(iElec)};
        
        if isempty(dat)
            continue;
        end
        
        %Reject bad trials
        good_trials_A = setdiff(find(dat.cond==1 & dat.trl_ind~=1),dat.bad_trials);
        good_trials_B = setdiff(find(dat.cond==2 & dat.trl_ind~=1),dat.bad_trials);
        good_trials_C = setdiff(find(dat.cond==3 & dat.trl_ind~=1),dat.bad_trials);
        
        %Retrieve HG envelope
        power_A = dat.power(good_trials_A,:);
        power_B = dat.power(good_trials_B,:);
        power_C = dat.power(good_trials_C,:);
        
        HG = [power_A ; power_B ; power_C];
        
        %Smooth data for visualization purposes
        for iTrial = 1:size(HG,1)
            HG(iTrial,:) = band_pass(HG(iTrial,:),srate,0,7);
        end
        
        
        h=axes('Position',gridPos{elec_locs{iGrid}(iElec)});
        %pcolor(tm_st:tm_en,1:size(HG,1),HG);
        pcolor(1:(bl+ps)*srate,1:size(HG,1),HG(:,1:(bl+ps)*srate));
        %  xlabel('ms'); ylabel('Z-Score');
        %  set(gca, 'XTick', [tm_st:jm:tm_en], 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');
        xlim([tm_st tm_en])
        % set(fgrid, 'PaperPositionMode', 'auto')
        hold on; plot([bl*srate bl*srate],[0 size(HG,1)],'-k','LineWidth',1)
        shading flat;
        caxis([-2 2]);
        
        % Plot electrode in title
        %title(num2str(elec))
        axis off;
        
        text(round((tm_en+tm_st)/2)-50,size(HG,1)+4,['e' num2str(es(iElec))],'FontSize',6)
        
        %make axis the electrode label, turn yaxis off
        %set(gca, 'XTick',round((tm_en+tm_st)/2),'XTickLabel', num2str(elec),'FontSize',5);
        %set(gca, 'YTickLabel',[]);
        cnt = cnt + 1;
    end
    
    
    
    saveas(gcf,[fig_path 'power_' grid_labels{iGrid} '.png'],'png')
    close
    
end
end


