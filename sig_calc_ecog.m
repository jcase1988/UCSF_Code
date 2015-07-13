clear

subj = 'EC70';
%blocks = {'B11', 'B12', 'B13'};
%blocks = {'B31' 'B33' 'B34'};
blocks = {'B11', 'B12', 'B13' 'B31' 'B33' 'B34'};
day = 'Day12';
srate = 400;
path = '/Users/johncase/Documents/UCSF Data/';
ana_path = [path subj '/' day '/analysis/'];
if ~exist(ana_path)
    mkdir(ana_path);
end
    
%% Define sliding window bins
bin_size_ms = 20; %bin size in ms
dt_ms = 10;       %sliding window increments in ms

bin_size = round(bin_size_ms/1000*srate);
dt = round(dt_ms/1000*srate);

bl_s = 0.5; %baseline in seconds
ps_s = 1;   %post-stim in seconds

bin_begs = 1:dt:round(((bl_s+ps_s)*srate)-bin_size); %from baseline to post-stim minus one window
bin_ends = bin_begs + bin_size;

[a stim_beg] = min(abs(bin_begs - (bl_s*srate)));

%onsets = round(ecogCAR.allstimtimes(:,1)*srate);
%offsets = round(ecogCAR.allstimtimes(:,2)*srate);

%% Load data from each block

for iBlock = 1:length(blocks)
    load([path subj '/' subj blocks{iBlock} '/analysis/HG_cond_power_z-scored.mat']) %HG_cond_power = {elecs}.(A/B/C) where A/B/C are the three different sounds
    trial_power_data{iBlock} = HG_cond_power;
    block_elecs{iBlock} = ~cellfun(@isempty,trial_power_data{iBlock});
    if iBlock == 1
        elecs = find(block_elecs{1});
    else
        elecs = intersect(elecs,find(block_elecs{iBlock}));
    end
end


%% Perform ANOVA tests on windows
sig_table = NaN(elecs(end),length(bin_begs));

for iElec = 1:length(elecs)
    elec = elecs(iElec);
    
    HG_A = []; HG_B = []; HG_C = [];
    for iBlock = 1:length(blocks)
        
        A_trials = trial_power_data{iBlock}{elec}.A;
        A_trials(find(trial_power_data{iBlock}{elec}.Ai==1),:) = []; %remove first trial
        HG_A = [HG_A ; A_trials];
        
        B_trials = trial_power_data{iBlock}{elec}.B;
        B_trials(find(trial_power_data{iBlock}{elec}.Bi==1),:) = []; %remove first trial
        HG_B = [HG_B ; B_trials];
        
        C_trials = trial_power_data{iBlock}{elec}.C;
        C_trials(find(trial_power_data{iBlock}{elec}.Ci==1),:) = []; %remove first trial
        HG_C = [HG_C ; C_trials];
        
    end
    
    for iBin = 1:length(bin_begs)
        clear p
        bin_beg = bin_begs(iBin);
        bin_end = bin_ends(iBin);
        
        A_bined = mean(HG_A(:,bin_beg:bin_end),2);
        B_bined = mean(HG_B(:,bin_beg:bin_end),2);
        C_bined = mean(HG_C(:,bin_beg:bin_end),2);
        
        ABC = [A_bined ; B_bined ; C_bined];
        group = [ones(length(A_bined),1) ; ones(length(B_bined),1)*2 ; ones(length(C_bined),1)*3];
        
        [p table stats] = anova1(ABC,group,'off');
        sig_table(elec,iBin) = p;
        
    end
    
end


%% count significant bins per elec

sig_counts = zeros(size(sig_table,1),1);
%[dummy stim_beg] = min(abs((bin_begs/srate)-bl_s)); %find first non-baseline bin

eoi = [0.1 0.75];
[d eoi_bin_beg] = min(abs((bin_begs/srate)-bl_s-eoi(1)));
[d eoi_bin_end] = min(abs((bin_begs/srate)-bl_s-eoi(2)));

for iElec = 1:size(sig_table,1)
    %  elec = elecs(iElec);
    sig_counts(iElec,1) =  length(find(sig_table(iElec,stim_beg:length(bin_begs)) < 0.05));
    
    % Find longest consectu
    x = sig_table(iElec,stim_beg:length(bin_begs)) < 0.05;
    sig_counts_consecutive(iElec,1) = max( diff( [0 (find( x== 0) )  numel(x) + 1] ) - 1);
end

save([ana_path 'sig_table_3_stimuli.mat'],'sig_table','sig_counts','sig_counts_consecutive','bin_begs','bin_size_ms','dt_ms');