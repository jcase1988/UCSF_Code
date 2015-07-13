clear

subj = 'EC77'; %subj = 'EC70';
%blocks = {'B11' 'B12' 'B13'};
%blocks = {'B31' 'B33' 'B34'};
%blocks = {'B11' 'B12' 'B13' 'B31' 'B33' 'B34'};'
blocks = {'B26' 'B27' 'B28' 'B42' 'B43' 'B44'};
day = 'Day1';
do_cons_offsets = 0; %if 1, time-lock to consonsant; otherwise, time-lock to stim

path = '/Users/johncase/Documents/UCSF Data/';

%load([path subj '/' day '/analysis/sig_table_3_stimuli.mat']);

%good_elecs = find(sig_counts>11)';
%good_elecs = find(sig_counts_consecutive>4)';


good_elecs = [243:244 227 211 209 195:196 178:179 149];

grid_elecs = 1:256;
good_elecs = good_elecs(ismember(good_elecs,grid_elecs));

srate = 400;

%onsets = round(ecogCAR.allstimtimes(:,1)*srate);
%offsets = round(ecogCAR.allstimtimes(:,2)*srate);

bl = 0.5; %baseline

%consonant offsets for the three stimuli
cons_offsets = [round(0.445*srate) round(0.319*srate) round(0.269*srate)];

%how long should the trial window be calculated?
if do_cons_offsets
    ps = 0.5;
    lock = 'Consonant';
else
    ps = 1;
    lock = 'Stimulus';
end

% These bin configs are used for the "dat_matrix" classification of <trials x electrodes*bins>
bin_size_2D = 0.050 * srate;
dt_2D = 0.010 * srate; %x in seconds time bins
bins_2D = round(1:dt_2D:((bl+ps)*srate)-bin_size_2D);

dat_matrix = []; % trials x (bins * electrodes)
dat_labels = {};

% These bin configs are used for the "dat_matrix_3D" classiizefication of <trials x electrodes x bins>
bin_size_3D = 0.050 * srate;
dt_3D = 0.010 * srate; %x in seconds time bins
bins_3D = round(1:dt_3D:((bl+ps)*srate)-bin_size_3D);

%plot x-tick at every 10th bin
plot_array = -bl*1000:(dt_3D/srate)*1000*10:ps*1000-((dt_3D/srate)*1000);
bin_samples = round(1000*((bins_3D/srate)-0.5));
for iTick = 1:length(plot_array)
    [a i] = min(abs(bin_samples-plot_array(iTick)));
    x_lab_cor(iTick) = i;
    plot_str{iTick} = num2str(plot_array(iTick));
end

%% Load z-scored data

for iBlock = 1:length(blocks)
    load([path subj '/' subj blocks{iBlock} '/analysis/HG_cond_power_z-scored.mat'])
    trial_power_data{iBlock} = HG_cond_power;
end

%% Arrange data into matrices

elec_cnt = 1;
elec_labels_2D = {}; elec_labels_3D = {}; dat_matrix_3D = []; dat_matrix = [];
for elec = good_elecs
    hgsignal = []; stim_labels = []; block_labels = [];
    for iBlock = 1:length(blocks)
        clear HG_cond_power
        dat = trial_power_data{iBlock}{elec};        
        
        %Reject bad trials
        good_trials_A = find(dat.cond==1 & dat.trl_ind~=1);
        good_trials_B = find(dat.cond==2 & dat.trl_ind~=1);
        good_trials_C = find(dat.cond==3 & dat.trl_ind~=1);
        
        %Retrieve HG envelope
        power_A = dat.power(good_trials_A,1:1.5*srate);
        power_B = dat.power(good_trials_B,1:1.5*srate);
        power_C = dat.power(good_trials_C,1:1.5*srate);        
        
        hgsignal = [hgsignal ; power_A ; power_B ; power_C];
        stim_labels = [stim_labels ; ones(length(good_trials_A),1)     ;
            ones(length(good_trials_B),1) * 2 ;
            ones(length(good_trials_C),1) * 3];
        
        block_labels = [block_labels ones(1,length(stim_labels))*iBlock];
        
    end
    
    
    %how much each trial binning should be offset
    %to time lock to consonants
    if do_cons_offsets
        wind_offsets = cons_offsets(stim_labels);
    else
        wind_offsets = zeros(1,length(stim_labels));
    end
    
    %Concatenate 2D matrix
    clear hg_bin
    for iBin = 1:length(bins_2D)
        
        hgsignal_temp = [];
        for iTrial = 1:length(stim_labels)
            wind = (bins_2D(iBin):bins_2D(iBin)+bin_size_2D-1)+wind_offsets(iTrial);
            hgsignal_temp = [hgsignal_temp ; hgsignal(iTrial,wind)];
        end
        
        hg_bin(:,iBin) = mean(hgsignal_temp,2);
        
        
        elec_labels_2D = [elec_labels_2D ['e' num2str(elec) '_' num2str(bins_2D(iBin))]];
    end
    dat_matrix = [dat_matrix hg_bin];
    
    %Concatenate 3D matrix
    clear hg_bin
    for iBin = 1:length(bins_3D)
        
        hgsignal_temp = [];
        for iTrial = 1:length(stim_labels)
            wind = (bins_3D(iBin):bins_3D(iBin)+bin_size_3D-1)+wind_offsets(iTrial);
            hgsignal_temp = [hgsignal_temp ; hgsignal(iTrial,wind)];
        end
        
        dat_matrix_3D(:,elec_cnt,iBin) = mean(hgsignal_temp,2);
        
        elec_labels_3D = [elec_labels_3D ['e' num2str(elec) '_' num2str(bins_3D(iBin))]];
    end
    
    elec_cnt = elec_cnt + 1;
end

dat_matrix_no_pac = dat_matrix;

%% K-means

idx = kmeans(dat_matrix_no_pac,3);
kmean_accuracy = length(find(idx == stim_labels))/length(stim_labels);


%% leave-one-out validation
clear accuracy

% Collapse
%%% stim_labels(find(stim_labels ~= 1)) = 2;
%%% stim_labels(find(stim_labels == 3)) = 2;


[dat_matrix score latent] = pca(dat_matrix_no_pac);
%dat_matrix = dat_matrix(:,1:5);

for i = 1:length(stim_labels)
    
    train_ind = setdiff(1:size(dat_matrix,1),i);
    test_ind = i;
    
    training = dat_matrix(train_ind,:);
    test = dat_matrix(test_ind,:);
    training_labels = stim_labels(train_ind)';
    test_labels = stim_labels(test_ind)';
    
    class = classify(test,training,training_labels);
    accuracy(i) = class==test_labels;
end

figure;
hist(accuracy)
title(['Mean = ' num2str(mean(accuracy)) ', SD = ' num2str(std(accuracy))])

%% classification per time bin
clear accuracy_time_course
for iBin = 1:length(bins_3D)
    clear accuracy
    for i = 1:length(stim_labels)
        
        %leave-one-out
        train_ind = setdiff(1:size(dat_matrix,1),i);
        test_ind = i;
        
        training = dat_matrix_3D(train_ind,:,iBin);
        test = dat_matrix_3D(test_ind,:,iBin);
        training_labels = stim_labels(train_ind)';
        test_labels = stim_labels(test_ind)';
        
        [class,err,POSTERIOR,logp,coeff] = classify(test,training,training_labels);
        accuracy(i) = class==test_labels;
    end
    
    accuracy_time_course(iBin) = mean(accuracy);
    
end

figure; plot(accuracy_time_course)
%set(gca, 'XTick', [1 51 100 150], 'XTickLabel', [-500 0 500 1000], 'XTickMode', 'manual', 'Layer', 'top');
xlim([1 length(accuracy_time_course)])
xlabel('MS')
ylabel('Accuracy')
title([lock '-locked Decoding']);
set(gca, 'XTick', x_lab_cor, 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');    

%% Investigate bin 106 for best features

% training = squeeze(dat_matrix_3D(:,:,106));
% test = squeeze(dat_matrix_3D(1,:,106)); %doesn't matter
% training_labels = stim_labels;
%
% [class,err,POSTERIOR,logp,coeff] = classify(test,training,training_labels);
%
% [a i12] = sort(abs(coeff(1,2).linear),'descend');
% [a i13] = sort(abs(coeff(1,3).linear),'descend');
% [a i23] = sort(abs(coeff(2,3).linear),'descend');
%
% be12 = good_elecs(i12(1:10));
% be13 = good_elecs(i13(1:10));
% be23 = good_elecs(i23(1:10));