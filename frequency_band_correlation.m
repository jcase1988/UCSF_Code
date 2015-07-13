clear
subj = 'EC82';
block = 'B49';

% Help find phase-providing central frequency and bandwidth that gives the
% best CFC (measured by MI) when coupled to high gamma (70-150 hz) across
% entire block

get_subj_globals(subj,block)
load([dtdir subj '_' block '_CAR.mat'])

freqs = 1:0.1:20;
bandwidth = 0.1;

elec = 310;

corr_mat = eye(length(freqs));

for f1 = 1:length(freqs)-1
    f1
    data1 = eegfilt(ecogCAR.data(elec,:),srate,freqs(f1),[]);
    data1 = eegfilt(data1,srate,[],freqs(f1)+bandwidth);
    
    for f2 = f1+1:length(freqs)
        data2 = eegfilt(ecogCAR.data(elec,:),srate,freqs(f2),[]);
        data2 = eegfilt(data2,srate,[],freqs(f2)+bandwidth);
        
        c = corrcoef(data1,data2);
        corr_mat(f1,f2) = c(1,2);
        corr_mat(f2,f1) = c(1,2);
    end
end