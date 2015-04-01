clear

subj = 'EC70';
blocks = {'B31' 'B34' 'B35'};
path = '/Users/johncase/Documents/UCSF Data/EC70';
dir = 'EC70B31';
CAR_datafile = 'EC70_B31_CAR.mat';
load([path '/' dir '/data/' CAR_datafile])

good_elecs = setdiff(ecogCAR.banks(1):ecogCAR.banks(end),ecogCAR.badChannels);
grid_elecs = 1:256;
good_elecs = good_elecs(ismember(good_elecs,grid_elecs));

good_elecs = [5, 6, 7, 22, 36, 38, 39, 40,  55, 56, 71, 96, 103, 112, 122, 134, 139, 149, 154, 165, 166, 181, 182, 196, 197, 229];

srate = ecogCAR.sampFreq;

onsets = round(ecogCAR.allstimtimes(:,1)*srate);
offsets = round(ecogCAR.allstimtimes(:,2)*srate);

bl = 0.5; %baseline
ps = 3;   %post-stim

dt = 0.050 * srate; %20 ms time bins
bins = round(((bl+0.25)*srate):dt:((bl+0.75)*srate)-dt);

dat_matrix = []; % trials x (bins * electrodes)
dat_labels = {};
cnt = 1;
for elec = good_elecs   
    clear erp band hg_signal hgsignal_z
    band = ecogCAR.data(elec,:);

        
    band = abs(my_hilbert(band, srate, 70, 150)).^2;
    

    ERP_A = []; ERP_B = []; ERP_C = [];
    for iTrial = 1:length(onsets)
        bl_onset = round(onsets(iTrial)-(bl*srate));
        bl_win = bl_onset:onsets(iTrial);
        win = bl_onset:round(onsets(iTrial)+(ps*srate));        
        hgsignal(iTrial,:) = band(win) - mean(band(bl_win));
        
        %z-transform
        bl_m = mean(band(bl_win));
        bl_sd = std(band(bl_win));
        zed = (band(win)-bl_m)./bl_sd;
        hgsignal_z(iTrial,:) = zed;
        
    end
    hgsignal(1,:) = [];
    hgsignal_z(1,:) = [];
    
    for iBin = 1:length(bins)
        hg_bin(:,iBin) = mean(hgsignal(:,bins(iBin):bins(iBin)+dt),2);
        dat_labels = [dat_labels ['e' num2str(elec) '_' num2str(bins(iBin))]];
    end
    
    dat_matrix = [dat_matrix hg_bin];
    
    erp_mean(:,cnt) = mean(hgsignal(:,(bl+0.25)*srate:(bl+0.75)*srate),2);
    erp_z_mean(:,cnt) = mean(hgsignal_z(:,(bl+0.25)*srate:(bl+0.75)*srate),2);
    cnt = cnt + 1;
end
labels = ecogCAR.stimID;

labels  = labels(2:end);


%% K-means

idx = kmeans(dat_matrix',2);
idx = idx';

%% LDA


% Make sure 50% of each category is in training set and test sets
clear accuracy
for i = 1:10000
    
    ind_a = shuffle(find(labels == 1));
    ind_b = shuffle(find(labels == 2));
    ind_c = shuffle(find(labels == 3));
        
    train_ind = [ind_a(1:round(length(ind_a)/2)) ...
                ind_b(1:round(length(ind_b)/2)) ...
                ind_c(1:round(length(ind_c)/2))];
            
    test_ind = setdiff(1:59,train_ind);
    
    training = erp_z_mean(train_ind,:);
    test = erp_z_mean(test_ind,:);
    training_labels = labels(train_ind)';
    test_labels = labels(test_ind)';
    
    class = classify(test,training,training_labels);
    accuracy(i) = length(find([class-test_labels]==0))/length(class);
end
figure; 
hist(accuracy)
title(['Mean = ' num2str(mean(accuracy)) ', SD = ' num2str(std(accuracy))])


%% leave-one-out validation
clear accuracy
for i = 1:length(labels)
    
    train_ind = setdiff(1:59,i);
    test_ind = i;
    
    training = erp_z_mean(train_ind,:);
    test = erp_z_mean(test_ind,:);
    training_labels = labels(train_ind)';
    test_labels = labels(test_ind)';
    
    class = classify(test,training,training_labels);
    accuracy(i) = class==test_labels;
end

figure; 
hist(accuracy)
title(['Mean = ' num2str(mean(accuracy)) ', SD = ' num2str(std(accuracy))])