clear
subj = 'EC77';
blocks = {'B26','B27','B28','B42','B43','B44'};

depth_banks = [257 276 ;
        277 308];

labels = {'AOF','POF','ITG','PTemp','AST','MST','PST','insula','sACC','iACC','amygdala','hippocampus'};
lab_ind = [257 261 265 271 277 281 285 289 293 297 301 305];
    
for iBlock = 1:length(blocks)
    clear dat_z dat_car dat_car_no_bank coh_z coh_car
    block = blocks{iBlock};
    get_subj_globals(subj,block);
    
    load([dtdir subj '_' block '.mat'])
    fig_path = [figdir 'raw_coherence/'];
    if ~exist(fig_path)
        mkdir(fig_path)
    end
    
    dat = ecogDS.data(257:308,:);
    dat = notchfilter(dat,400,60); dat = notchfilter(dat,400,120); dat = notchfilter(dat,400,180);
    
    m_dat = mean(dat,2);
    sd_dat = std(dat,0,2);
    
    for i = 1:length(dat)
        dat_z(:,i) = (dat(:,i)-m_dat) ./ sd_dat;
    end
    
    %No bank CAR
    CAR = mean(dat_z(setdiff(1:size(dat_z,1),badChannels-256),:),1);
    for i = 1:length(dat)
        dat_car_no_bank(:,i) = dat_z(:,i) - CAR(i);
    end
    
    %Bank CAR
    clear CAR
    for iBank = 1:size(depth_banks,1)
        elecs = (depth_banks(iBank,1)-256):(depth_banks(iBank,2)-256);
        good_elecs = setdiff(elecs,badChannels-256);
        
        CAR = mean(dat_z(good_elecs,:),1);
        for i = 1:length(dat)
            dat_car(elecs,i) = dat_z(elecs,i) - CAR(i);
        end
    end
    
    
    x = hilbert(dat_z);
    phases = angle(x);
    
    for elec1 = 1:size(dat_z,1)
        for elec2 = 1:size(dat_z,1)
            coh_z(elec1,elec2) = abs(sum(exp(1i * (phases(elec1,:) - phases(elec2,:))), 'double'))...
                / size(phases,2);
        end
    end
    

    
    x = hilbert(dat_car_no_bank);
    phases = angle(x);
    
    for elec1 = 1:size(dat_car_no_bank,1)
        for elec2 = 1:size(dat_car_no_bank,1)
            coh_car_no_bank(elec1,elec2) = abs(sum(exp(1i * (phases(elec1,:) - phases(elec2,:))), 'double'))...
                / size(phases,2);
        end
    end
    
    
    x = hilbert(dat_car);
    phases = angle(x);
    
    for elec1 = 1:size(dat_car,1)
        for elec2 = 1:size(dat_car,1)
            coh_car(elec1,elec2) = abs(sum(exp(1i * (phases(elec1,:) - phases(elec2,:))), 'double'))...
                / size(phases,2);
        end
    end
    
    figure;
    subplot(3,2,1); heatmap(coh_z); caxis([0 0.5]); title('coherence Raw'); colorbar
    set(gca,'YTick',lab_ind-256,'YTickLabel',labels)
    subplot(3,2,2);  hist(reshape(coh_z,1,size(coh_z,1)*size(coh_z,2)),100); title('coherence Raw');% xlim([0 0.3])    
    
    %figure;
    subplot(3,2,3);  heatmap(coh_car_no_bank); caxis([0 0.5]); title('coherence CAR No Bank'); colorbar    
    set(gca,'YTick',lab_ind-256,'YTickLabel',labels)
    subplot(3,2,4);  hist(reshape(coh_car_no_bank,1,size(coh_car_no_bank,1)*size(coh_car_no_bank,2)),100); title('coherence CAR No Bank');% xlim([0 0.3])
    
    %figure;
    subplot(3,2,5);  heatmap(coh_car); caxis([0 0.5]); title('coherence CAR'); colorbar    
    set(gca,'YTick',lab_ind-256,'YTickLabel',labels)
    subplot(3,2,6);  hist(reshape(coh_car,1,size(coh_car,1)*size(coh_car,2)),100); title('coherence CAR');% xlim([0 0.3])
    
    saveas(gcf,[fig_path 'CAR_all_depths.fig']);
end