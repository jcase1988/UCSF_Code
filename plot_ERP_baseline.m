function plot_ERP_baseline(subj,block,power)

path = ['/Users/johncase/Documents/UCSF Data/' subj '/' subj block];
CAR_path = [path '/data/' subj '_' block '_CAR.mat'];
load(CAR_path)

fig_path = [path '/figures/'];
ana_path = [path '/analysis/'];

if power
    fig_path = [fig_path 'CARed_HG_power/power_'];
else
    fig_path = [fig_path 'CARed_ERP/erp_'];
end

if ~exist(fig_path)
    mkdir(fig_path)
end

good_elecs = setdiff(ecogCAR.banks(1):ecogCAR.banks(end),ecogCAR.badChannels);
srate = ecogCAR.sampFreq;


onsets = round(ecogCAR.allstimtimes(:,1)*srate);
offsets = round(ecogCAR.allstimtimes(:,2)*srate);

bl = 0.5; %baseline
ps = 1;   %post-stim

start_time_window = -bl*1000;
tm_st = round( start_time_window ./1000 *srate);
tm_en = round( ps*srate);
plot_jump = 250;
jm = round(plot_jump./1000*srate);

for elec = good_elecs
    clear erp band
    band = ecogCAR.data(elec,:);
    
    
    if power
        band = abs(my_hilbert(band, srate, 70, 150)).^2;
    end
    
    ERP_A = []; ERP_B = []; ERP_C = []; A_ind = []; B_ind = []; C_ind = [];
    for iTrial = 1:length(onsets)
        bl_onset = round(onsets(iTrial)-(bl*srate));
        bl_win = bl_onset:onsets(iTrial);
        win = bl_onset:round(onsets(iTrial)+(ps*srate));
        erp(iTrial,:) = band(win) - mean(band(bl_win));
        
        %z-transform
        bl_m = mean(band(bl_win));
        bl_sd = std(band(bl_win));
        zed = (band(win)-bl_m)./bl_sd;
        erp_z(iTrial,:) = zed;
        
        if ecogCAR.stimID(iTrial) == 1
            ERP_A = [ERP_A ; zed];
            A_ind = [A_ind ; iTrial]; %include index so that I can exclude first trial later
        elseif ecogCAR.stimID(iTrial) == 2
            ERP_B = [ERP_B ; zed];
            B_ind = [B_ind ; iTrial];
        elseif ecogCAR.stimID(iTrial) == 3
            ERP_C = [ERP_C ; zed];
            C_ind = [C_ind ; iTrial];
        end
    end
    
    HG_cond_power{elec}.A = ERP_A;
    HG_cond_power{elec}.B = ERP_B;
    HG_cond_power{elec}.C = ERP_C;
    HG_cond_power{elec}.Ai = A_ind;
    HG_cond_power{elec}.Bi = B_ind;
    HG_cond_power{elec}.Ci = C_ind;
    
    mean_erp = mean(erp,1);
    mean_erp_z = mean(erp_z,1);
    ERP_A = mean(ERP_A,1);
    ERP_B = mean(ERP_B,1);
    ERP_C = mean(ERP_C,1);
    
    for z = 1:length([tm_st:jm:tm_en])
        plot_str{z} = round(start_time_window+(z-1)*plot_jump); 
    end
    
    % Plot signal averaged over all conditions
    fgrid = figure('visible','off');
    plot(tm_st:tm_en,band_pass(mean_erp,srate,0,7),'LineWidth',3); 
    xlabel('ms'); ylabel('uV');
    set(gca, 'XTick', [tm_st:jm:tm_en], 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');
    xlim([tm_st tm_en])    
    set(fgrid, 'PaperPositionMode', 'auto')
    saveas(gcf, [fig_path 'e' num2str(elec) '.jpg'], 'jpg')
    close
    
    % Plot Z-scored signal averaged over all conditions
    figure('visible','off');
    plot(band_pass(mean_erp_z,srate,0,7)); 
    xlabel('ms'); ylabel('Z-Score');
    set(gca, 'XTick', [tm_st:jm:tm_en], 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');
    xlim([tm_st tm_en])    
    set(fgrid, 'PaperPositionMode', 'auto')
    saveas(gcf, [fig_path 'Z-score_e' num2str(elec) '.jpg'], 'jpg')
    close
    
    % Plot Z-scored signal for each condition
    fgrid = figure('visible','off');
    plot(tm_st:tm_en,band_pass(ERP_A,srate,0,7),'LineWidth', 2)
    hold on;
    plot(tm_st:tm_en,band_pass(ERP_B,srate,0,7),'-r','LineWidth', 2)
    hold on;
    plot(tm_st:tm_en,band_pass(ERP_C,srate,0,7),'-g','LineWidth', 2)
    set(gca, 'XTick', [tm_st:jm:tm_en], 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');
    legend('aagaa','iyfiy','uwshuw')
    xlabel('ms'); 
    ylabel('Z-Score'); 
    xlim([tm_st tm_en])    
    set(fgrid, 'PaperPositionMode', 'auto')
    
    saveas(gcf, [fig_path 'Z-score_conds_e' num2str(elec) '.jpg'], 'jpg')
    close;
end

save([ana_path 'HG_cond_power_z-scored.mat'],'HG_cond_power');