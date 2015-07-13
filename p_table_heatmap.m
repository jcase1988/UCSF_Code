function p_table_heatmap(p_table,subj,block,time_lock,beg_win,en_win,dt,bin_size,fig_label)
%time_lock = 0 for stimulus lock
%            1 for consonant lock             


%%

elec_locs = [16:16:256 ; 15:16:255 ; 14:16:254 ; 13:16:253 ; 12:16:252 ; 11:16:251 ; 10:16:250 ; 9:16:249 ; 8:16:248 ; 7:16:247 ; 6:16:246 ; 5:16:245 ; 4:16:244 ; 3:16:243 ; 2:16:242 ; 1:16:241];

if strcmp('EC82',subj)
    elecs_locs_depth = [257:266 ; 267:276 ; 277:280 0 0 0 0 0 0 ; 281:284 0 0 0 0 0 0 ; 285:288 0 0 0 0 0 0 ; 289:292 0 0 0 0 0 0 ; 293:296 0 0 0 0 0 0 ; 297:302 0 0 0 0 ; 303:306 0 0 0 0 0 0 ; 307:310 0 0 0 0 0 0 ; 311:314 0 0 0 0 0 0];
    elecs_locs_depth_lab = {'FG', 'FG', 'Amyg', 'Hippo', 'AST', 'MST', 'PST', 'ITG','F1','sOFC','iOFC'};
elseif strcmp('EC77',subj)
    elecs_locs_depth = [257:260 0 0 ; 261:264 0 0; 265:270 ; 271:276 ; 277:280 0 0; 281:284 0 0; 285:288 0 0; 289:292 0 0; 293:296 0 0; 297:300 0 0; 301:304 0 0; 305:308 0 0];
    elecs_locs_depth_lab = {'AOF','POF','ITG','Post-temp','AST','MST','PST','Insula','sACC','iACC','Amy','Hippo'};
    
end


if strcmp(subj,'EC77') %if right
    elec_locs = fliplr(elec_locs);
end

get_subj_globals(subj,block{1})

%for figure file name
bins = num2cell((beg_win:dt:en_win-bin_size)*1000);

fig_path = [subjdir [block{:}] '/figures/anova_heatmaps/'];
if ~exist(fig_path)
    mkdir(fig_path)
end

ext = {'' 'stim-ME', 'value-ME','stim-value-X'};


if time_lock == 1
    lock = 'consonant';
elseif time_lock == 0
    lock = 'stimulus';
end

for iComp = 2:4 %2: stim main-effect, 3: value main-effect, 4: stim/value interaction
    
    fig_path_ext = [fig_path fig_label '/' ext{iComp} '/'];
    if ~exist(fig_path_ext)
        mkdir(fig_path_ext)
    end
    
    %parse all of the p-value data
    for iBin = 1:length(p_table)        
        %lateral grid
        all_dat(:,iBin) = p_table{iBin}(iComp,:);
    end
    
    %delete all p-value data that isn't < 0.01 for three consecutive bins
    n_consec = 3;
    p_cutoff = 0.01;
    thes_dat = ones(size(all_dat));
    for iElec = 1:size(all_dat,1)
        for iBin = 1:size(all_dat,2)-1
            
            if all_dat(iElec,iBin)<p_cutoff
                cnt = 1;
                for iBin2 = iBin+1:size(all_dat,2)
                    if all_dat(iElec,iBin2)<p_cutoff
                        cnt = cnt + 1;
                        continue;
                    else
                        if cnt >= n_consec
                            thes_dat(iElec,iBin:iBin2-1) = all_dat(iElec,iBin:iBin2-1);
                        end
                        break;
                    end
                end
            end
        end
    end
                    
    for iBin = 1:length(p_table) 
        dat = thes_dat(:,iBin);
        
        figure('visible','off');
        heatmap(dat(elec_locs),[],[],elec_locs);
        caxis([0 0.012]);
        colorbar;       
        
        title([num2str(bins{iBin}) ' to ' num2str(bins{iBin}+bin_size*1000)]);
        
        saveas(gcf,[fig_path_ext lock '_' (num2str(round(bins{iBin}))) 'ms.png'],'png');
        close;
        
        %depths
        for iRow = 1:size(elecs_locs_depth,1)
            for iCol = 1:size(elecs_locs_depth,2)
                elec = elecs_locs_depth(iRow,iCol);
                if elec == 0 || ismember(elec,badChannels)
                    depth_dat(iRow,iCol) = 1;
                else
                    depth_dat(iRow,iCol) = thes_dat(elec,iBin);
                end
            end
        end
        figure('visible','off');
        heatmap(depth_dat,[],elecs_locs_depth_lab,elecs_locs_depth);
        caxis([0 0.012]);
        colorbar;   
        
        title([num2str(bins{iBin}) ' to ' num2str(bins{iBin}+bin_size*1000)]);
        
        saveas(gcf,[fig_path_ext 'depths_' lock '_' (num2str(round(bins{iBin}))) 'ms.png'],'png');
        close;
        
    end
end

