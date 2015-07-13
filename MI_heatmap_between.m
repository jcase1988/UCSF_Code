clear
subj = 'EC71';
block = {'B2','B11'};

fp = 1:0.1:15.1;
fp_bandwidth = 0.5:0.1:5.1;

phase_elec = 1:60;
amp_elec =   90;
%84 85 100 162

if iscell(block)
    
    get_subj_globals(subj,block{1})
    figure_path = [subjdir [block{:}] '/figures/MI_scan/'];
    if ~exist(figure_path)
        mkdir(figure_path)
    end

    for iElec_amp = 1:length(amp_elec)
        for iElec_phase = 1:length(phase_elec)
       
        clear MI_all
        h =figure('Position', [1 1 2400 1260]);
        for iBlock = 1:length(block)
            clear MI
            get_subj_globals(subj,block{iBlock})

            clear MI
            load([anadir 'MI/e' num2str(phase_elec(iElec_phase)) '_e' num2str(amp_elec(iElec_amp)) '.mat'])            
            MI_all(:,:,iBlock) = MI;
            
        end

        %plot 
        %c = [0 max(max(max(MI_all)))];        
        
        for iBlock = 1:length(block)
            subplot(1,length(block)+1,iBlock);
            heatmap(squeeze(MI_all(:,:,iBlock)),fp_bandwidth,fp);            
            caxis([0 0.0001]);
            title(['e' num2str(phase_elec(iElec_phase)) '  e' num2str(amp_elec(iElec_amp)) ' ' block{iBlock}])
            colorbar;
        end
        
        subplot(1,length(block)+1,length(block)+1);
        heatmap(squeeze(mean(MI_all,3)),fp_bandwidth,fp);
    %    caxis([0 c(2)]);
        
        colorbar;
        title('Mean')
           
        saveas(h,[figure_path 'between_e'  num2str(phase_elec(iElec_phase)) '_e' num2str(amp_elec(iElec_amp)) '.fig'],'fig')
        saveas(h,[figure_path 'between_e'  num2str(phase_elec(iElec_phase)) '_e' num2str(amp_elec(iElec_amp)) '.png'],'png')
        close
        end
    end
else


    for iElec_amp = 1:length(amp_elec)
        for iElec_phase = 1:length(phase_elec)

            get_subj_globals(subj,block)

            figure_path = [figdir 'MI_scan/'];
            if ~exist(figure_path)
                mkdir(figure_path)
            end

            clear MI
            load([anadir '/MI/e' num2str(phase_elec(iElec_phase)) '_e' num2str(amp_elec(iElec_amp)) '.mat'])

            MI_all(:,:,iElec_phase) = MI;
            
            %h = figure('visible','off');
            h =figure;
            heatmap(MI,fp_bandwidth,fp);
            title(['e' num2str(phase_elec(iElec_phase)) '  e' num2str(amp_elec(iElec_amp))])
            colorbar; 
            saveas(h,[figure_path 'between_e'  num2str(phase_elec(iElec_phase)) '_e' num2str(amp_elec(iElec_amp)) '.fig'],'fig')
            saveas(h,[figure_path 'between_e'  num2str(phase_elec(iElec_phase)) '_e' num2str(amp_elec(iElec_amp)) '.png'],'png')
            close(h)
        end
    end
end