clear
subj = 'EC71';
block = 'B2';
get_subj_globals(subj,block);

anin_to_use = 2; 
threshold = 0.3;
wpath = [UCSF_dir '/ecog preprocessing/Erin/RLsounds'];
load([wpath '/names.mat'])
tasktype = 2;

[eventtimes,eventnames,aga,ushu,ifi,clicks]=DetectEvents_Erin(dtdir,wpath,anin_to_use,names,tasktype,threshold);
'Done'


%% For all blocks

allstimtimes = eventtimes;

stimID(strcmp(eventnames,'aagaa')) = 1;
stimID(strcmp(eventnames,'iyfiy')) = 2;
stimID(strcmp(eventnames,'uwshuw')) = 3;
stimID(strcmp(eventnames,'click')) = 10;    

add_subj_globals(subj,block,stimID,allstimtimes);
%% For learning blocks

feedbacktimes = [];
sample1times = [];
sample2times = [];
clicktimes = [];

cnt = 1;
for i = 1:4:size(eventtimes,1)
    sample1times = [sample1times ; eventtimes(i,:)];
    sample2times = [sample2times ; eventtimes(i+1,:)];
    clicktimes = [clicktimes ; eventtimes(i+2,:)];
    feedbacktimes = [feedbacktimes ; eventtimes(i+3,:)];
    
    offerval(cnt,1:2) = [stimID(i) stimID(i+1)];
    feedbackval(cnt) = stimID(i+3);
    cnt = cnt + 1;
end
    
add_subj_globals(subj,block,feedbacktimes,sample1times,sample2times,clicktimes,offerval,feedbackval);