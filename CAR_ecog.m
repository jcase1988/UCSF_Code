clear
subj = 'EC70';
block = 'B34';

path = ['/Users/johncase/Documents/UCSF Data/' subj '/' subj block];
CAR_path = [path '/data/' subj '_' block '_CAR.mat'];
raw_path = [path '/data/' subj '_' block '.mat'];
load(raw_path)


if exist(CAR_path)
    overwrite = input('CAR file already exists. Overwrite (1) or quit(2)?');
    if overwrite ~= 1
        return
    end
end

%notch filters
ecogDS.data = notchfilter(ecogDS.data,ecogDS.sampFreq,60);
ecogDS.data = notchfilter(ecogDS.data,ecogDS.sampFreq,120);
ecogDS.data = notchfilter(ecogDS.data,ecogDS.sampFreq,180);

banks = default_banks();
    
%Perform a bank-wise CAR
ecogCAR = ecogDS;
ecogCAR.banks = banks;
ecogCAR.data = zeros(banks(length(banks),2),size(ecogDS.data,2));
for iBank = 1:size(banks,1)
    clear CAR
    elecs = banks(iBank,1):banks(iBank,2);
    good_elecs = setdiff(elecs,ecogDS.badChannels); %exclude bad channels from CAR and referencing
    CAR_elecs = setdiff(good_elecs,ecogDS.spikeChannels); %exclude spike channels from CAR only
    CAR = sum(ecogDS.data(CAR_elecs,:),1)/length(CAR_elecs);
    
    if length(CAR_elecs) == 1
        ecogCAR.data(CAR_elecs,:) = ecogDS.data(CAR_elecs,:);
        good_elecs = setdiff(good_elecs,CAR_elecs); %if there is only one good elec in bank, do not re-reference that channel        
    end
    
    for e = 1:length(good_elecs)
        ecogCAR.data(good_elecs(e),:) = ecogDS.data(good_elecs(e),:) - CAR;
    end
end



save([path '/' dir '/data/' datafile '_CAR.mat'],'ecogCAR')


