%function scan_freq_for_MI(subj,block)
clear
subj = 'EC71';
blocks = {'B2' 'B11'};

for iBlock = 1:length(blocks)

    figure;
    
    block = blocks{iBlock};
    
    % Contrast MI of three stimuli during baseline block

    get_subj_globals(subj,block)
    load([dtdir subj '_' block '_CAR.mat'])
    cols = {'-b','-r','-g'};

    fps = [2.9 6.5 4.7 4.3 5.5 5.1 5.2 4 4.7];
    fp_bandwidths = [2.7 4.9 2.7 4.9 5 2.2 4.4 1.7 4.6];

    fps = ones(1,100)*6;
    fp_bandwidths = ones(1,100)*4;
    
    
    %offset_all = [0:0.1:2]; 
    %ps = 0.5; %time window = onsets of stimulus to "ps" secs later

    offset_all = 0;
    ps = 3;

    amp_elec = 90;
    phase_elecs = [1:48];

    for iElec = 1:length(phase_elecs)

        phase_elec = phase_elecs(iElec);
        
        %phase-providing frequency
        fp = fps(iElec);
        fp_bandwidth = fp_bandwidths(iElec);

        %amplitude-providing frequency
        fa = [70 150];

        %define phase bins
        n_bins = 100;
        bin_size = 2*pi/n_bins;
        bins = -pi:bin_size:(pi-bin_size);


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
        % pow = pow(setdiff(t_win,bad_samp));

        %Calculate MI for each phase-providing central-frequency / bandwidth


        %Do phase-providing phase filtering
        phase = zeros(1,length(ecogCAR.data));
        phase = eegfilt(ecogCAR.data(phase_elec,:), srate, [], fp+(fp_bandwidth/2));
        phase = eegfilt(phase, srate, fp-(fp_bandwidth/2), []);
        phase = angle(hilbert(phase));
        %   phase = phase(setdiff(t_win,bad_samp));

        
        for iOffset = 1:length(offset_all)

            offset = offset_all(iOffset);

            for iStim = 1:3
                onsets = round(allstimtimes(find(stimID==iStim),1)*srate);
                ind = [];
                for onset = onsets'
                    ind = [ind onset+round(offset*srate):onset+round((offset+ps)*srate)];
                end

                ind = setdiff(ind,bad_samp);



                %Calculate mean amplitude within each phase bin to yield a
                %distribution of amplitude(phase)
                bin_dist = zeros(1,length(bins));
                for iBin = 1:length(bins)

                    bin_ind = find(phase>=bins(iBin) & phase<bins(iBin)+bin_size);
                    bin_ind = bin_ind(find(ismember(bin_ind,ind)));
                    bin_dist(iBin) = mean(pow(bin_ind));
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
                MI(iElec,iBlock,iStim) = (log(length(bins)) - h_p) / log(length(bins));
                hold on; plot(bin_dist,cols{iOffset});
            end
        end
    end
end