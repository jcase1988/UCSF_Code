clear
subj = 'EC70';
block = 'B13';

dirpath = ['/Users/johncase/Documents/UCSF Data/' subj '/' subj block '/data/'];
load([dirpath subj '_' block '_CAR.mat'])
srate = ecogCAR.sampFreq;

%create events structure
for iTrial = 1:length(ecogCAR.allstimtimes)
    events(iTrial).type='111';
    events(iTrial).latency=round(ecogCAR.allstimtimes(iTrial)*srate);
    events(iTrial).duration = round(ecogCAR.allstimtimes(iTrial,2)*srate) - round(ecogCAR.allstimtimes(iTrial,1)*srate);
end

for iElec = 1:size(ecogCAR.data,1)
    HG(iElec,:) = abs(my_hilbert(ecogCAR.data(iElec,:), srate, 70, 150)).^2;
end

eegplot(HG,'srate',ecogCAR.sampFreq,'events',events)

if ~isfield(ecogCAR,'bad_epochs')
    bad_epochs = cell(size(ecogCAR.data,1),1);
else
    bad_epochs = ecogCAR.bad_channel_epochs;
end

inp = '';
all_inps = {};
while true
    inp = input('','s');
    
    if strcmp(inp,'q')
        break;
    end
    
    all_inps = [all_inps inp];
    
end

for iInp = 1:length(all_inps)
    inp = all_inps{iInp};
    inp_split = strsplit(inp,' ');
    
    beg_epoch = str2num(inp_split{1});
    end_epoch = str2num(inp_split{2});
    
    for elec = elecs
        bad_epochs{elec} = [bad_epochs{elec} beg_epoch end_epoch];
    end
end

ecogCAR.bad_epochs = bad_epochs;

sav = input('Save bad epoch? Yes (1) or No (2)');
if sav == 1
    save([dirpath subj '_' block '_CAR.mat'],'ecogCAR')
end
    