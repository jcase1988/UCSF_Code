clear
subj = 'EC70';
block = 'B11';

anin_to_use = 1; 
dtpath = ['/Users/johncase/Documents/UCSF Data/' subj '/' subj block '/'];
wpath = '/Users/johncase/Documents/UCSF Data/ecog preprocessing/Erin/RLsounds';
load([wpath '/names.mat'])
tasktype = 2;

[eventtimes,eventnames,aga,ushu,ifi,clicks]=DetectEvents_Erin(dtpath,wpath,anin_to_use,names,tasktype);