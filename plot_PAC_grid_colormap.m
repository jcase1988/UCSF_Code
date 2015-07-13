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


fgrid = figure;
set(fgrid, 'Position', ScSz);
for elec = elecs
    if isempty(PAC_cond{elec_locs(elec)})
        continue;
    end
    
    subplot(16,16,elec);
    for iCond = 1:3
        inds = setdiff(find(PAC_cond{elec_locs(elec)}.cond==iCond),   PAC_cond{elec_locs(elec)}.bad_trials);
        PAC_dat = PAC_cond{elec_locs(elec)}.power;
        wind = round((bl+st)*srate):round((bl+en)*srate);
        
        dat(iCond) = mean(abs(mean(PAC_dat(inds,wind),1)));
    end
    
    norm_multiple = 1/max(dat);
    
    bar(1,1,'FaceColor',dat*norm_multiple)
    %bar(1,1,'FaceColor',dat)
    
      %  bar(iCond,dat,col{iCond});
      %  hold on;
      %  ylim([0 0.6])
        set(gca,'XTickLabel','','YTickLabel','')
      %  title(num2str(elec_locs(elec)));
        
    end
end
