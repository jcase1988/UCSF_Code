clear
subj = 'EC71';
block = 'B2';

fp = 1:0.1:15.1;
fp_bandwidth = 0.5:0.1:5.1;

%phase_elec = [18 19 49 50 51 66 67 82 83 84 85 100 101 102 115 116 117 131 132 133 146 147 148 149 162 163 164 165 179 180 181 211 277 278 279 280 281 282 283 284 307 308 309 310 311 312 313 314];
%amp_elec =   [18 19 49 50 51 66 67 82 83 84 85 100 101 102 115 116 117 131 132 133 146 147 148 149 162 163 164 165 179 180 181 211 277 278 279 280 281 282 283 284 307 308 309 310 311 312 313 314];
%phase_elec = [18 19 49 50 51 66 67 82 83 84 85 100 101 102 115 116 117 131];
%amp_elec =   [18 19 49 50 51 66 67 82 83 84 85 100 101 102 115 116 117 131];

phase_elec = 1:64;
amp_elec = 90;


if iscell(block)
    
    get_subj_globals(subj,block{1})
    figure_path = [subjdir [block{:}] '/figures/MI_scan/'];
    if ~exist(figure_path)
        mkdir(figure_path)
    end

    for iElec = 1:length(amp_elec)  
        clear MI_all
        h =figure('Position', [1 1 2400 1260]);
        for iBlock = 1:length(block)
            clear MI
            get_subj_globals(subj,block{iBlock})

            clear MI
            load([anadir 'MI/e' num2str(phase_elec(iElec)) '_e' num2str(amp_elec(iElec)) '.mat'])            
            MI_all(:,:,iBlock) = MI;
            
        end

        %plot 
        c = [0 max(max(max(MI_all(10:end,:,:))))];        
        
        for iBlock = 1:length(block)
            subplot(1,length(block)+1,iBlock);
            heatmap(squeeze(MI_all(:,:,iBlock)),fp_bandwidth,fp);
            caxis(c);
            title(['e' num2str(amp_elec(iElec)) ' ' block{iBlock}])
            colorbar;
        end
        
        subplot(1,length(block)+1,length(block)+1);
        heatmap(squeeze(mean(MI_all,3)),fp_bandwidth,fp);
        caxis([0 c(2)]);
        
        colorbar;
        title('Mean')
           
        saveas(h,[figure_path 'within_e' num2str(phase_elec(iElec)) '.fig'],'fig')
        close
    end
else


    for iElec = 1:length(amp_elec)

        get_subj_globals(subj,block)

        figure_path = [figdir 'MI_scan/'];
        if ~exist(figure_path)
            mkdir(figure_path)
        end


        clear MI
        load([anadir 'MI/e' num2str(phase_elec(iElec)) '_e' num2str(amp_elec(iElec)) '.mat'])

        MI_all_elecs(:,:,amp_elec(iElec)) = MI;
        
        %h = figure('visible','off');
        h =figure;
        heatmap(MI,fp_bandwidth,fp);
        title(['e' num2str(phase_elec(iElec)) ' e' num2str(amp_elec(iElec))])
        colorbar; 
        saveas(h,[figure_path 'within_e' num2str(phase_elec(iElec)) '.fig'],'fig')
        close(h)
    end
end