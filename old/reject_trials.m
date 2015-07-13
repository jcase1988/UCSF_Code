clear
subj = 'EC70';
block = 'B11';

CAR_path = ['/Users/johncase/Documents/UCSF Data/' subj '/' subj block '/data/' subj '_' block '_CAR.mat'];
load(CAR_path)

if ~isfield(ecogCAR,'bad_epochs')
    ecogCAR.bad_epochs = [];
end

while true
    
    temp = pop_eegplot(ecogCAR.data,'srate',ecogCAR.sampFreq,'dispchans',64);
    uiwait; 
    
    chans = input('Which channels do these bad epochs correspond to? Input as vector or "0" to apply to all channels.\n');
    if chans == 0
        chans = size(ecogCAR.data,1);
    end
    
    % Channels x bad epoch index x (begin/end samples)
    for e = chans
        ecogCAR.bad_epochs(e,:,1:2) = [squeeze(ecogCAR.bad_epochs(e,:,1:2)) TMPREJ(:,1:2)];
    end
    
    cond = input('Input 1 to reject more trials.\n');
    if ~(cond == 1)
        break;
    end
    
end

save([path '/' dir '/data/' CAR_datafile],'ecogCAR')
