function [seg_bled] = baseline_trial(seg,srate,bl)       
    
    bl_win = round(1:(bl*srate)-1);
    
    
    for e = 1:size(seg,1)
        for iTrial = 1:size(seg,3)
            
            seg_bled(e,:,iTrial) = seg(e,:,iTrial)-mean(seg(e,bl_win,iTrial));
        end
    end
end
            