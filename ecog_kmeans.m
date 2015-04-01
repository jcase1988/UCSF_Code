clear

subj = 'EC70';
blocks = {'B31' 'B33' 'B34'};

path = '/Users/johncase/Documents/UCSF Data/';


%good_elecs = setdiff(ecogCAR.banks(1):ecogCAR.banks(end),ecogCAR.badChannels);
%grid_elecs = 1:256;
%good_elecs = good_elecs(ismember(good_elecs,grid_elecs));


%good_elecs = [6 7 50 54 96 128 134 135];
%good_elecs = [6 7];
good_elecs = [6 7 19 20 50 53 54 57 60 96 114 119 128 134 135 187 195 197 207 256 259 264 278 282];

srate = 400;

%onsets = round(ecogCAR.allstimtimes(:,1)*srate);
%offsets = round(ecogCAR.allstimtimes(:,2)*srate);

bl = 0.5; %baseline
ps = 1;   %post-stim

% These bin configs are used for the "dat_matrix" classification of <trials x electrodes*bins>
bin_size_2D = 0.050 * srate; 
dt_2D = 0.010 * srate; %x in seconds time bins
bins_2D = round(((bl+0.1)*srate):dt_2D:((bl+0.4)*srate)-bin_size_2D);

dat_matrix = []; % trials x (bins * electrodes)
dat_labels = {};


% These bin configs are used for the "dat_matrix_3D" classiizefication of <trials x electrodes x bins>
bin_size_3D = 0.050 * srate; 
dt_3D = 0.010 * srate; %x in seconds time bins
bins_3D = round(1:dt_3D:((bl+1)*srate)-bin_size_3D);


%% 

for iBlock = 1:length(blocks)
    load([path subj '/' subj blocks{iBlock} '/analysis/HG_cond_power.mat'])
    trial_power_data{iBlock} = HG_cond_power;
end

%%
elec_cnt = 1;
elec_labels_2D = {}; elec_labels_3D = {}; dat_matrix_3D = []; dat_matrix = [];
for elec = good_elecs  
    hgsignal = []; stim_labels = []; trial_num_labels = [];
    for iBlock = 1:length(blocks)
        clear HG_cond_power
        HG_cond_power = trial_power_data{iBlock};
        hgsignal = [hgsignal ; HG_cond_power{elec}.A ; HG_cond_power{elec}.B ; HG_cond_power{elec}.C];
        stim_labels = [stim_labels ; ones(size(HG_cond_power{elec}.A,1),1)     ; 
                           ones(size(HG_cond_power{elec}.B,1),1) * 2 ; 
                           ones(size(HG_cond_power{elec}.C,1),1) * 3];
                       
        trial_num_labels = [trial_num_labels ; HG_cond_power{elec}.Ai     ;
            HG_cond_power{elec}.Bi ;
            HG_cond_power{elec}.Ci ];
 
    end
        
    %remove the first trial of each block becuase auditory response is
    %likely to be an outlier
    remove_first = find(trial_num_labels == 1);
    hgsignal(remove_first,:) = [];
    stim_labels(remove_first) = [];
    
    %Concatenate 2D matrix
    clear hg_bin        
    for iBin = 1:length(bins_2D)
        hg_bin(:,iBin) = mean(hgsignal(:,bins_2D(iBin):bins_2D(iBin)+bin_size_2D-1),2);
        elec_labels_2D = [elec_labels_2D ['e' num2str(elec) '_' num2str(bins_2D(iBin))]];
    end    
    dat_matrix = [dat_matrix hg_bin];           
    
    %Concatenate 3D matrix
    clear hg_bin        
    for iBin = 1:length(bins_3D)
        dat_matrix_3D(:,elec_cnt,iBin) = mean(hgsignal(:,bins_3D(iBin):bins_3D(iBin)+bin_size_3D-1),2);
        elec_labels_3D = [elec_labels_3D ['e' num2str(elec) '_' num2str(bins_3D(iBin))]];
    end    
    
    elec_cnt = elec_cnt + 1;
end


%% K-means

idx = kmeans(dat_matrix,3);
kmean_accuracy = length(find(idx == stim_labels))/length(stim_labels);

%% LDA
% 
% 
% % Make sure 50% of each category is in training set and test sets
% clear accuracy
% for i = 1:10000
%     
%     ind_a = shuffle(find(stim_labels == 1));
%     ind_b = shuffle(find(stim_labels == 2));
%     ind_c = shuffle(find(stim_labels == 3));
%         
%     train_ind = [ind_a(1:round(length(ind_a)/2)) ...
%                 ind_b(1:round(length(ind_b)/2)) ...
%                 ind_c(1:round(length(ind_c)/2))];
%             
%     test_ind = setdiff(1:size(dat_matrix,1),train_ind);
%     
%     training = erp_z_mean(train_ind,:);
%     test = erp_z_mean(test_ind,:);
%     training_labels = labels(train_ind)';
%     test_labels = labels(test_ind)';
%     
%     class = classify(test,training,training_labels);
%     accuracy(i) = length(find([class-test_labels]==0))/length(class);
% end
% figure; 
% hist(accuracy)
% title(['Mean = ' num2str(mean(accuracy)) ', SD = ' num2str(std(accuracy))])


%% leave-one-out validation
clear accuracy

 stim_labels(find(stim_labels ~= 3)) = 1;
 stim_labels(find(stim_labels == 3)) = 2;

 dat_matrix_old = dat_matrix;
 [dat_matrix score latent] = pca(dat_matrix);
 dat_matrix = dat_matrix(:,1:3);

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

% figure; 
% hist(accuracy)
% title(['Mean = ' num2str(mean(accuracy)) ', SD = ' num2str(std(accuracy))])

%% classification per time bin
clear accuracy_time_course
for iBin = 1:length(bins_3D)
    clear accuracy
    for i = 1:length(stim_labels)
        
        train_ind = setdiff(1:size(dat_matrix,1),i);
        test_ind = i;
        
        training = dat_matrix_3D(train_ind,:,iBin);
        test = dat_matrix_3D(test_ind,:,iBin);
        training_labels = stim_labels(train_ind)';
        test_labels = stim_labels(test_ind)';
        
        class = classify(test,training,training_labels);
        accuracy(i) = class==test_labels;
    end
    
    accuracy_time_course(iBin) = mean(accuracy);
    
end

figure; plot(accuracy_time_course)
set(gca, 'XTick', [1 51 100 150], 'XTickLabel', [-500 0 500 1000], 'XTickMode', 'manual', 'Layer', 'top');