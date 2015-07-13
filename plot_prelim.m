clear

subj = 'EC82';
block = 'B47';

get_subj_globals(subj,block);
path = ['/Users/johncase/Documents/UCSF Data/' subj '/' subj block];
raw_path = [path '/data/' subj '_' block '.mat'];
CAR_path = [path '/data/' subj '_' block '_CAR.mat'];
load(raw_path)
load(CAR_path)

fig_path = [path '/figures/'];


%banks = ecogCAR.banks;
    
%Plot spectrograms of raw and CARed signals
if ~exist([fig_path 'raw_spec/'])
    mkdir(fig_path)
    mkdir([fig_path 'raw_spec/'])
    mkdir([fig_path 'CARed_spec/'])
    mkdir([fig_path 'raw_CARed_per/'])
    mkdir([fig_path 'raw_CARed_ERP/'])
    
end

seg = segment_trial(ecogCAR,0.5,2,allstimtimes);
erps = squeeze(mean(seg,3));
good_elecs = setdiff(banks(1):banks(end),ecogDS.badChannels);
clear per_raw per_CARed
for e = 1:length(good_elecs)
    [per_raw(:,good_elecs(e)) x f] = periodogram(ecogDS.data(good_elecs(e),:),[],200,400);
    [per_CARed(:,good_elecs(e)) x f] = periodogram(ecogCAR.data(good_elecs(e),:),[],200,400);
end

for e = 1:length(good_elecs)
     % Plot raw spectragram
     figure('visible','off'); 
     specgram(ecogDS.data(good_elecs(e),:),[],400); caxis([1 100])
     saveas(gcf, [fig_path 'raw_spec/spec_e' num2str(good_elecs(e)) '.jpg'], 'jpg')
     close
     
     % Plot CARed spectragram 
      figure('visible','off'); 
      specgram(ecogCAR.data(good_elecs(e),:),[],400); caxis([1 100])
      saveas(gcf, [fig_path 'CARed_spec/spec_e' num2str(good_elecs(e)) '.jpg'], 'jpg')
      close
    
     % Plot periodograme
    figure('visible','off'); 
    semilogy(x,per_raw(:,good_elecs(e)))
    hold on;
    semilogy(x,per_CARed(:,good_elecs(e)),'-r')
    ylim([0 100000])
    legend('Raw','CARed')
    saveas(gcf, [fig_path 'raw_CARed_per/per_e' num2str(good_elecs(e)) '.jpg'], 'jpg')
    close
    
   % Plot CAR ERP
   figure('visible','off'); 
   plot(erps(good_elecs(e),:))
   saveas(gcf, [fig_path 'raw_CARed_ERP/per_e' num2str(good_elecs(e)) '.jpg'], 'jpg')
    close
end