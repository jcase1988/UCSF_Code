elec = 101;
%%

dat = HG_cond_power{elec};
%values = setdiff(unique(dat.value),-999)';
cons_offsets = [round(0.445*srate) round(0.319*srate) round(0.269*srate)];
bl = 0.5 * srate;
labels = {'aagaa','iyfiy','uwshuw'};

st_ms = -0.6;
en_ms = 0.4;
step_ms = 0.1;

st_sa = round(st_ms*srate);
en_sa = round(en_ms*srate);
step_sa = round(step_ms*srate);

%for 
x_labels = num2cell(st_ms*1000:100:en_ms*1000);
x_ticks = 1:0.1*srate:(en_sa-st_sa);

for iVal = 1:length(values)
    tm_win = (bl+cons_offsets(iVal)+st_sa):(bl+cons_offsets(iVal)+en_sa);
    
    inds = find(dat.value == values(iVal))';
    val_pow = dat.power(inds,:);
    avg_ERP = mean(dat.power(inds,tm_win),1);
    
    
    for i = 1:size(val_pow)
        [c l] = xcorr(val_pow(i,tm_win),avg_ERP,40);
        [a max_corr_ind] = max(c);
        lags(i) = l(max_corr_ind);
    end
    
    clear fitted_pow
    for iTrial = 1:size(val_pow)
        fitted_pow(iTrial,:) = val_pow(iTrial,tm_win+lags(iTrial));
    end
    
    fitted_pow = mean(fitted_pow,1);
    
    %band-pass
  %  avg_ERP = band_pass(avg_ERP,400,0,20);
  %  fitted_pow = band_pass(fitted_pow,400,0,20);
    
    figure; plot(avg_ERP)
    hold on; plot(fitted_pow,'-r')
    legend('HG','HG Fitted')
    ylabel('Z-score')
    xlabel('Ms from Consonant')
    xlim([0 length(avg_filt)])
    title(['electrode ' num2str(elec) ' - ' labels{iVal}])
    set(gca,'XTick',x_ticks,'XTickLabel',x_labels);
end