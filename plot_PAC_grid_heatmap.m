function plot_PAC_grid(subj,block)
% get_subj_globals(subj,block)
% 
% % Reassign paths if "All" or "Day" Variable
% if strcmp(block(1:3),'All')
%     anadir = anadir_all;
%     figdir = figdir_all;
% elseif strcmp(block(1:3),'Day')
%     anadir = anadir_day;
%     figdir = figdir_day;
% end
% 
% fig_path = [figdir 'power_grid_aud_resp/'];
% 
% if ~exist(fig_path)
%     mkdir(fig_path)
% end
% 
% good_elecs = setdiff(banks(1):banks(end),badChannels);

ScSz=[1 1 2400 1260];


col = {'b','r','g'};
stim_label = {'aagaa','iyfiy','uwshuw'};

bl = 0.5;

%begin and end time window relative to simulus onset
st = 0;
en = 0.5;

if iscell(block)
    for iBlock = 1:length(block)
        get_subj_globals(subj,block{iBlock})
        load([anadir 'PAC_cond_theta_gamma.mat'])
        PAC{iBlock} = PAC_cond;
    end
    
    for elec = 1:length(PAC{1})
        PAC_cat{elec}.power = [];
        PAC_cat{elec}.cond = [];
        PAC_cat{elec}.value = [];
        PAC_cat{elec}.trl_ind = [];
        PAC_cat{elec}.bad_trials = [];
        
        trl_cnt = 0;
        for iBlock = 1:length(block)
            if ~isempty(PAC{iBlock}{elec})
                PAC_cat{elec}.power = [PAC_cat{elec}.power ; PAC{iBlock}{elec}.power];
                PAC_cat{elec}.cond = [PAC_cat{elec}.cond ; PAC{iBlock}{elec}.cond];
                PAC_cat{elec}.value = [PAC_cat{elec}.value ; PAC{iBlock}{elec}.value];
                PAC_cat{elec}.trl_ind = [PAC_cat{elec}.trl_ind ; PAC{iBlock}{elec}.trl_ind];
                
                PAC_cat{elec}.bad_trials = [PAC_cat{elec}.bad_trials  PAC{iBlock}{elec}.bad_trials+trl_cnt];
                trl_cnt = trl_cnt + size(PAC_cat{elec}.power,1);
            end
        end
        if isempty(PAC_cat{elec}.power)
            PAC_cat(elec) = [];
        end
    end
        
    PAC_cond = PAC_cat;
        
else
    load([anadir 'PAC_cond_theta_gamma.mat'])
end

elecs = 1:256;
elec_locs = [1:16 ; 17:32 ; 33:48 ; 49:64 ; 65:80 ; 81:96 ; 97:112 ; 113:128 ; 129:144 ; 145:160 ; 161:176 ; 177:192 ; 193:208 ; 209:224 ; 225:240 ; 241:256];
elec_locs = flipud(flipud(elec_locs)')';
elec_locs = reshape(elec_locs,1,256);

for iCond = 1:3

    fgrid = figure;
    set(fgrid, 'Position', ScSz);
    for elec = elecs
        if isempty(PAC_cond{elec_locs(elec)})
            dat(elec) = 0;
            continue;
        end

        %subplot(16,16,elec);

            inds = setdiff(find(PAC_cond{elec_locs(elec)}.cond==iCond),   PAC_cond{elec_locs(elec)}.bad_trials);
            PAC_dat = PAC_cond{elec_locs(elec)}.power;
            wind = round((bl+st)*srate):round((bl+en)*srate);

            dat(elec) = mean(abs(mean(PAC_dat(inds,wind),1)));
    end

    cnt = 1;
   for irow = 1:16
       for icol = 1:16
         dat_mat(irow,icol) = dat(cnt);
         cnt = cnt + 1;
       end
   end
    
   heatmap(dat_mat)      
    
    caxis([0 0.4])
   
   title(stim_label{iCond});
    
    
end

