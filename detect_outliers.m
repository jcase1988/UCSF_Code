clear
subj = 'EC77';
block = 'B26';

path = ['/Users/johncase/Documents/UCSF Data/' subj '/' subj block];
CAR_path = [path '/data/' subj '_' block '_CAR.mat'];
ana_path = [path '/analysis/'];
load(CAR_path)
load([ana_path '/HG_cond_power_z-scored.mat'])

good_elecs = setdiff(ecogCAR.banks(1):ecogCAR.banks(end),ecogCAR.badChannels);

trial_counts = zeros(size(ecogCAR.allstimtimes,1),1);
elec_counts = zeros(size(ecogCAR.data,1),1);
for elec = good_elecs
    dat = HG_cond_power{elec};
    A_m = mean(dat.A(:,200:end),2);
    B_m = mean(dat.B(:,200:end),2);
    C_m = mean(dat.C(:,200:end),2);
    
    A_std = std(A_m);
    B_std = std(B_m);
    C_std = std(C_m);
    
    A_ind = dat.Ai;
    B_ind = dat.Bi;
    C_ind = dat.Ci;
    
    outliers{elec} = [A_ind(find(mean(A_m)>A_m+A_std*3))
                      A_ind(find(mean(A_m)<A_m-A_std*3))
                
                      B_ind(find(mean(B_m)>B_m+B_std*3))
                      B_ind(find(mean(B_m)<B_m-B_std*3))
                
                      C_ind(find(mean(C_m)>C_m+C_std*3))
                      C_ind(find(mean(C_m)<C_m-C_std*3))];
    for trl = outliers{elec}
        trial_counts(trl) = trial_counts(trl) + 1;
    end
    elec_counts(elec) = length(outliers{elec});
end