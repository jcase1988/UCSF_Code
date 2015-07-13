%function scan_freq_for_MI(subj,block)
clear
subj = 'EC71';
block = 'B11';

% Help find phase-providing central frequency and bandwidth that gives the
% best CFC (measured by MI) when coupled to high gamma (70-150 hz) across
% entire block

get_subj_globals(subj,block)
load([dtdir subj '_' block '_CAR.mat'])

phase_elec = 39;

for amp_elec = 39

    %phase-providing frequency
    %fp = 1:0.1:15;            %try 1-20 hz in steps of 0.1
    %fp_bandwidth = 0.5:0.1:5; %adjust bandwidth

    fp = 5.7;
    fp_bandwidth = 4.6;
    
    %amplitude-providing frequency
    fa = [70 150]; 

    %define phase bins
    n_bins = 20;
    bin_size = 2*pi/n_bins;
    bins = -pi:bin_size:(pi-bin_size);

 
    %define time_window (roughly entire block) (exclude artifacts samps later)
    t_0 = round(allstimtimes(1,1)*srate);
    t_end = round((allstimtimes(end,1)+4) * srate);
    t_win = t_0:t_end;


    %identify artifact samples
    bad_samp = [];
    for iEpoch = 1:size(per_chan_bad_epochs{phase_elec},1)
        beg = per_chan_bad_epochs{phase_elec}(iEpoch,1)*srate;
        en = per_chan_bad_epochs{phase_elec}(iEpoch,2)*srate;

        bad_samp = [bad_samp beg:en];
    end

    if phase_elec ~= amp_elec
        for iEpoch = 1:size(per_chan_bad_epochs{amp_elec},1)
            beg = per_chan_bad_epochs{amp_elec}(iEpoch,1)*srate;
            en = per_chan_bad_epochs{amp_elec}(iEpoch,2)*srate;

            bad_samp = unique([bad_samp beg:en]);
        end
    end

    %Do high-gamma filtering
    pow = eegfilt(ecogCAR.data(amp_elec,:), srate, fa(1), []);
    pow = eegfilt(pow, srate, [], fa(2));
    pow = abs(hilbert(pow));
    pow = pow(setdiff(t_win,bad_samp));

    %Calculate MI for each phase-providing central-frequency / bandwidth
    MI = zeros(length(fp),length(fp_bandwidth));
    for iFreq = 1:length(fp)
        for iBand = 1:length(fp_bandwidth)

            if fp(iFreq)-(fp_bandwidth(iBand)/2) < 0.5 %if the lower bound of freq range is less than 0.5
                MI(iFreq,iBand) = 0;
                continue;
            end

            disp(fp(iFreq))
            disp(fp_bandwidth(iBand))

            %Do phase-providing phase filtering
            phase = zeros(1,length(ecogCAR.data));
            phase = eegfilt(ecogCAR.data(phase_elec,:), srate, [], fp(iFreq)+(fp_bandwidth(iBand)/2));
            phase = eegfilt(phase, srate, fp(iFreq)-(fp_bandwidth(iBand)/2), []);        
            phase = angle(hilbert(phase));
            phase = phase(setdiff(t_win,bad_samp));

            %Calculate mean amplitude within each phase bin to yield a
            %distribution of amplitude(phase)
            bin_dist = zeros(1,length(bins));
            for iBin = 1:length(bins)
                bin_dist(iBin) = mean(pow(phase>=bins(iBin) & phase<bins(iBin)+bin_size));
            end

            %Normalize distribution to yield pseudo "probability density function" (PDF)
            bin_dist = bin_dist ./ sum(bin_dist);

            %calculate Shannon entropy of PDF
            h_p = 0;
            for iBin = 1:length(bins)
                h_p = h_p - bin_dist(iBin)*log(bin_dist(iBin));
            end

            %MI = (Kullback-Leibler distance between h_p and uniform
            %distribution) / (Entropy of uniform distribution) (see
            %http://jn.physiology.org/content/104/2/1195)
            MI(iFreq,iBand) = (log(length(bins)) - h_p) / log(length(bins));


        end
    end

    if exist([anadir 'MI_scan.mat'])
        load([anadir 'MI_scan.mat'])
    end
    MI_scan.([phase_elec '_' amp_elec]) = MI;
    save([anadir 'MI_scan.mat'],'MI_scan')
end