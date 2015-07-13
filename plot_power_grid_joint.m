%function plot_power_grid(subj,block,power)
clear
subj = 'EC77';
%blocks = {'B11', 'B12', 'B13'};
%blocks = {'B31', 'B33', 'B34'};
%blocks = {'B11', 'B12', 'B13', 'B31', 'B33', 'B34'};
blocks = {'B26','B27','B28'};
day = 'Day1';
power = 1;
srate = 400;

path = ['/Users/johncase/Documents/UCSF Data/' subj '/' day '/'];
ana_path = [path 'analysis/'];
fig_path_grid = [path 'figures/single_trial_plots/'];
fig_path_line = [path 'figures/line_plots/'];


if ~exist(fig_path_grid)
    mkdir(fig_path_grid)
end
if ~exist(fig_path_line)
    mkdir(fig_path_line)
end
if ~exist(ana_path)
    mkdir(ana_path)
end


%% Define sliding window bins
bin_size_ms = 20; %bin size in ms
dt_ms = 10;       %sliding window increments in ms 

bin_size = round(bin_size_ms/1000*srate);
dt = round(dt_ms/1000*srate);

bl = 0.5; %baseline in seconds
ps = 1;   %post-stim in seconds

bin_begs = 1:dt:round(((bl+ps)*srate)-bin_size); %from baseline to post-stim minus one window
bin_ends = bin_begs + bin_size;

% What is the index of the closest bin to the onset of the stimulus?
[a stim_bin] = min(abs((bin_begs/400)-0.5));

%onsets = round(ecogCAR.allstimtimes(:,1)*srate);
%offsets = round(ecogCAR.allstimtimes(:,2)*srate);

if exist([ana_path 'HG_cond_power_z-scored_binned.mat'])
    load([ana_path 'HG_cond_power_z-scored_binned.mat']);
else
    
    HG_A = []; HG_B = []; HG_C = []; HG = [];
    onset_total = 0;
    onset_conds = [0 0 0];
    
    for iBlock = 1:length(blocks)
        z_score_path = ['/Users/johncase/Documents/UCSF Data/' subj '/' subj blocks{iBlock} '/analysis/HG_cond_power_z-scored.mat'];
        load(z_score_path)
        
        HG_data{iBlock} = HG_cond_power;
        
        good_elecs_blocks{iBlock} = find(~cellfun(@isempty,HG_cond_power,'UniformOutput',1));        
        
        if iBlock == 1
            good_elecs = good_elecs_blocks{1};
        else
            good_elecs = intersect(good_elecs,good_elecs_blocks{iBlock});
        end
        
        onset_total = onset_total + length(HG_cond_power{good_elecs_blocks{iBlock}(1)}.Ai) + length(HG_cond_power{good_elecs_blocks{iBlock}(1)}.Bi) + length(HG_cond_power{good_elecs_blocks{iBlock}(1)}.Ci);
        onset_conds = onset_conds + [length(HG_cond_power{good_elecs_blocks{iBlock}(1)}.Ai)  length(HG_cond_power{good_elecs_blocks{iBlock}(1)}.Bi)  length(HG_cond_power{good_elecs_blocks{iBlock}(1)}.Ci)];
        
    end
    


    %%


    start_time_window = -bl*1000;
    tm_st = round( start_time_window ./1000 *srate);
    tm_en = round( ps*srate);
    plot_jump = 250;
    jm = round(plot_jump./1000*srate);

    %initialize data matrices
    HG_all = zeros(good_elecs(end),onset_total,length(bin_begs));
    HG_A = zeros(good_elecs(end),onset_conds(1),length(bin_begs));
    HG_B = zeros(good_elecs(end),onset_conds(2),length(bin_begs));
    HG_C = zeros(good_elecs(end),onset_conds(3),length(bin_begs));

    for elec = good_elecs

        cnt = 1; cnt_A = 1; cnt_B = 1; cnt_C = 1;
        for iBlock = 1:length(blocks)
            clear erp erp_z band
            HG_unbinned_all = [HG_data{iBlock}{elec}.A ; HG_data{iBlock}{elec}.B ; HG_data{iBlock}{elec}.C];                                

            % All conditions
            for iTrial = 1:size(HG_unbinned_all,1)

                band = band_pass(HG_unbinned_all(iTrial,:),srate,0,7);
                % Average within time bins
                for iBin = 1:length(bin_begs)
                    zed_bined = mean(band(bin_begs(iBin):bin_ends(iBin)),2);
                    HG_all(elec,cnt,iBin) = zed_bined;
                end

                cnt = cnt + 1;            
            end     

            % A condition
            for iTrial = 1:size(HG_data{iBlock}{elec}.A,1)

                band = band_pass(HG_data{iBlock}{elec}.A(iTrial,:),srate,0,7);
                % Average within time bins
                for iBin = 1:length(bin_begs)
                    zed_bined = mean(band(bin_begs(iBin):bin_ends(iBin)),2);
                    HG_A(elec,cnt_A,iBin) = zed_bined;
                end

                cnt_A = cnt_A + 1;            
            end    

            % B condition
            for iTrial = 1:size(HG_data{iBlock}{elec}.B,1)

                band = band_pass(HG_data{iBlock}{elec}.B(iTrial,:),srate,0,7);
                % Average within time bins
                for iBin = 1:length(bin_begs)
                    zed_bined = mean(band(bin_begs(iBin):bin_ends(iBin)),2);
                    HG_B(elec,cnt_B,iBin) = zed_bined;
                end

                cnt_B = cnt_B + 1;            
            end  

            % C condition
            for iTrial = 1:size(HG_data{iBlock}{elec}.C,1)

                band = band_pass(HG_data{iBlock}{elec}.C(iTrial,:),srate,0,7);
                % Average within time bins
                for iBin = 1:length(bin_begs)
                    zed_bined = mean(band(bin_begs(iBin):bin_ends(iBin)),2);
                    HG_C(elec,cnt_C,iBin) = zed_bined;
                end

                cnt_C = cnt_C + 1;            
            end  


        end
    end
    save([ana_path 'HG_cond_power_z-scored_binned.mat'],'HG_all','HG_A','HG_B','HG_C');
end %end does not exist HG_cond_power_z-scored_binned

%% Plot grid

elecs = 1:256;
elec_locs = [1:16 ; 17:32 ; 33:48 ; 49:64 ; 65:80 ; 81:96 ; 97:112 ; 113:128 ; 129:144 ; 145:160 ; 161:176 ; 177:192 ; 193:208 ; 209:224 ; 225:240 ; 241:256];
elec_locs = flipud(flipud(elec_locs)')';
elec_locs = reshape(elec_locs,1,256);

ScSz=[1 1 2400 1260];
pos_grid=[1 1 2400 1260];

%ScSz=[1 1 1680 1050];
%pos_grid=[1 1 1680 1050];

gridPos = DivideScreen(16,16, ScSz, 20, 30, pos_grid); % Change first two variables to your first grid/strip's dimensions.


fgrid = figure('visible','off');
set(fgrid, 'Position', ScSz,'color', [1 1 1], 'MenuBar', 'figure');

for elec = elecs
    
    h=axes('Position',gridPos{elec_locs(elec)});
    pcolor(1:length(bin_begs),1:size(HG_all,2),squeeze(HG_all(elec,:,:)))
    %  xlabel('ms'); ylabel('Z-Score');
    %  set(gca, 'XTick', [tm_st:jm:tm_en], 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');
    xlim([1 length(bin_begs)])
    % set(fgrid, 'PaperPositionMode', 'auto')
    hold on; plot([stim_bin stim_bin],[0 size(HG_all,2)],'-k','LineWidth',1) %Plot stimulus line
    shading flat;
    caxis([-2 2]);
    
    % Plot electrode in title
    %title(num2str(elec))
    axis off;
    
    text(round(length(bin_begs)/2),size(HG_all,2)+14,['e' num2str(elec)],'FontSize',5)
    
    %make axis the electrode label, turn yaxis off
    %set(gca, 'XTick',round((tm_en+tm_st)/2),'XTickLabel', num2str(elec),'FontSize',5);
    %set(gca, 'YTickLabele(HG)',[]);
    
end

saveas(gcf,[fig_path_grid 'z-scored_grid.png'],'png')


%%

xticks = {'-400' '-200' '0' '200' '400' '600' '800' '1000'};
for iTick = 1:length(xticks)
    tick = str2num(xticks{iTick});
    [a q] = min(abs((bin_begs/srate)-bl-tick/1000));
    xtick_loc(iTick) = q;
end

for elec = good_elecs

% Plot Z-scored signal for each condition
    fgrid = figure('visible','off');
    plot(1:length(bin_begs),squeeze(mean(HG_A(elec,:,:),2)),'LineWidth', 2)
    hold on;
    plot(1:length(bin_begs),squeeze(mean(HG_B(elec,:,:),2)),'-r','LineWidth', 2)
    hold on;
    plot(1:length(bin_begs),squeeze(mean(HG_C(elec,:,:),2)),'-g','LineWidth', 2)
    set(gca, 'XTick', xtick_loc, 'XTickLabel', xticks, 'XTickMode', 'manual', 'Layer', 'top');
    legend('aagaa','iyfiy','uwshuw')
    xlabel('ms'); 
    ylabel('Z-Score'); 
    xlim([1 length(bin_begs)])    
    set(fgrid, 'PaperPositionMode', 'auto')
    
    saveas(gcf, [fig_path_line 'By_conditions_e' num2str(elec) '.jpg'], 'jpg')
    close;

    % Plot Z-scored signal for each condition
    fgrid = figure('visible','off');
    plot(1:length(bin_begs),squeeze(mean(HG_all(elec,:,:),2)),'LineWidth', 2)   
    set(gca, 'XTick', xtick_loc, 'XTickLabel', xticks, 'XTickMode', 'manual', 'Layer', 'top');    
    xlabel('ms'); 
    ylabel('Z-Score'); 
    xlim([1 length(bin_begs)])    
    set(fgrid, 'PaperPositionMode', 'auto')
    
    saveas(gcf, [fig_path_line 'Mean_over_trials_e' num2str(elec) '.jpg'], 'jpg')
    close;

end
