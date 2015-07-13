%% Second script to run after receiving data (after defining banks in default_banks.m)
clear
subj = 'EC71';
block = 'B2';
day = 'Day1';
task = 'RL_baseline'; %'RL_baseline' 'RL_learning'
srate = 400;
%% Input bad/spike channels and define banks

badChannels = [49 82 161:512];
spikeChannels = [141:144];
[banks bank_labels] = default_banks(subj);

%% Elec labels
% 
% bank_labels = {};
% bank_labels(1:256) = {'Grid'};
% bank_labels(257:276) = {'FG'};
% bank_labels(277:280) = {'Amyg'};
% bank_labels(281:284) = {'Hippo'};
% bank_labels(285:288) = {'ATP1'};
% bank_labels(289:292) = {'MST'};
% bank_labels(293:296) = {'PST'};
% bank_labels(297:302) = {'ITG'};
% bank_labels(303:306) = {'F1'};
% bank_labels(307:310) = {'sOF'};
% bank_labels(311:314) = {'iOF'};

%% Add globals

% Setup paths
add_subj_globals(subj,block) 

% Add additional information
add_subj_globals(subj,block,day);
add_subj_globals(subj,block,task,badChannels,spikeChannels,banks,srate);