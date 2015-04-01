%% Use this script to look determine if the banks are appropriate. If a small subset of electrodes
%% have a high variance, exclude them from your bank

clear
subj = 'EC70';
block = 'B33';
path = ['/Users/johncase/Documents/UCSF Data/' subj];
dir = [subj block];
datafile = [subj '_' block '.mat'];

load([path '/' dir '/data/' datafile])
srate = ecogDS.sampFreq;

%notch filters
ecogDS.data = notchfilter(ecogDS.data,ecogDS.sampFreq,60);
ecogDS.data = notchfilter(ecogDS.data,ecogDS.sampFreq,120);
ecogDS.data = notchfilter(ecogDS.data,ecogDS.sampFreq,180);

banks = default_banks();

bad_channels = [];
for iBank = 1:size(banks,1)
    eegplot(ecogDS.data(banks(iBank,1):banks(iBank,2),:),'srate',srate)
    uiwait;
    bad_channels = [bad_channels input('Bad channels to exclude?\n') + (banks(iBank,1)-1)];    
    elecs_exclude{iBank,1} = input('Additional electrodes to exclude from this bank?\n');   
    elecs_exclude{iBank,1} = elecs_exclude{iBank,1} + (banks(iBank,1)-1);
end