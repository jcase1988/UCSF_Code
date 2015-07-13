function [beta, p] = regression_statevals_helper(subj,ecogCAR,badChannels,bl_onsets,onsets,bl,ps,per_chan_bad_epochs,stateval,fig_path,titl)

    %% Setup plot parameters    
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

    %Load Data
    srate = ecogCAR.sampFreq;
    

    good_elecs = setdiff(1:size(ecogCAR.data,1),badChannels);
    p(badChannels) = 1;
    beta(badChannels) = 0;
    for elec = good_elecs
        signal = ecogCAR.data(elec,:);
        pow = [];
        sta = [];
        for iTrial = 1:length(onsets)
            bad_trial = 0;
            win = [bl_onsets(iTrial)-bl onsets(iTrial)+ps(2)];


            %identify bad trials
            for iEpoch = 1:size(per_chan_bad_epochs{elec},1)
                beg = per_chan_bad_epochs{elec}(iEpoch,1)*srate;
                en = per_chan_bad_epochs{elec}(iEpoch,2)*srate;

                if (beg > win(1) && beg < win(end)) || (en > win(1) && en < win(end))
                    bad_trial = 1;
                end
            end

            if bad_trial            
                continue;
            end

            bl_win = round((((bl_onsets(iTrial)-bl)*srate) : (bl_onsets(iTrial)*srate)));
            ps_win = round((((onsets(iTrial)+ps(1))*srate) : ((onsets(iTrial)+ps(2))*srate)));

            bl_m = mean(signal(bl_win));
            bl_sd = std(signal(bl_win));
            zed = (signal(ps_win)-bl_m)/bl_sd;

            pow = [pow mean(zed)];
            sta = [sta stateval(iTrial)];
        end

        [b dev stats] = glmfit(sta,pow);
        beta(elec) = b(2);
        p(elec) = stats.p(2);
    end
    
    
    %Plot
    
    figure('visible','off');
    heatmap(p(elec_locs),[],[],elec_locs);
    colorbar;
    caxis([0 0.08])
    title([titl ' p']);
    saveas(gcf,[fig_path '_grid_p.png'],'png');

    figure('visible','off');
    heatmap(beta(elec_locs),[],[],elec_locs);
    colorbar;   
    caxis([-0.1 0.1])
    title([titl ' Beta']);
    saveas(gcf,[fig_path '_grid_beta.png'],'png');
    close
    
    %depths
    for iRow = 1:size(elecs_locs_depth,1)
        for iCol = 1:size(elecs_locs_depth,2)
            elec = elecs_locs_depth(iRow,iCol);
            if elec == 0 || ismember(elec,badChannels)
                depth_p(iRow,iCol) = 1;
                depth_beta(iRow,iCol) = 0;
            else
                depth_p(iRow,iCol) = p(elec);
                depth_beta(iRow,iCol) = beta(elec);
            end
        end
    end
    figure('visible','off');
    heatmap(depth_p,[],elecs_locs_depth_lab,elecs_locs_depth);
    colorbar;
    caxis([0 0.08])
    title([titl ' p']);
    saveas(gcf,[fig_path '_depths_p.png'],'png');
    close
    
    figure('visible','off'); 
    heatmap(depth_beta,[],elecs_locs_depth_lab,elecs_locs_depth);
    colorbar;
    caxis([-0.1 0.1])
    title([titl ' Beta']);
    saveas(gcf,[fig_path '_depths_beta.png'],'png');
    close
    
    
end
