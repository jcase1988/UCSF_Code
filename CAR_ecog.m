function CAR_ecog(subj,block,varargin)

    get_subj_globals(subj,block);

    CAR_path = [dtdir subj '_' block '_CAR.mat'];
    raw_path = [dtdir subj '_' block '.mat'];
    load(raw_path)


    if exist(CAR_path) && isempty(varargin)
        overwrite = input('CAR file already exists. Overwrite (1) or quit(2)?\n');
        if overwrite ~= 1
            return
        end
    end

    %notch filters
    ecogDS = notf(ecogDS);

    if ~exist('banks')
        banks = default_banks(subj);
        add_subj_globals(subj,block,banks);
    end

    %Perform a bank-wise CAR
    ecogCAR = ecogDS;
    ecogCAR.data = zeros(banks(length(banks),2),size(ecogDS.data,2));
    for iBank = 1:size(banks,1)
        clear CAR
        elecs = banks(iBank,1):banks(iBank,2);
        
        %if depths, z-score signal
        if (strcmp(subj,'EC77') && elecs(1) > 256) || (strcmp(subj,'EC82') && elecs(1)>276) || (strcmp(subj,'EC71') && elecs(1)>128)
            for e = elecs
                m = mean(ecogDS.data(e,:));
                sd = std(ecogDS.data(e,:));
                ecogDS.data(e,:) = (ecogDS.data(e,:)-m)./sd;
            end
        end
        
        good_elecs = setdiff(elecs,badChannels); %exclude bad channels from CAR and referencing
        CAR_elecs = setdiff(good_elecs,spikeChannels); %exclude spike channels from CAR only
        CAR = sum(ecogDS.data(CAR_elecs,:),1)/length(CAR_elecs);

        if length(CAR_elecs) <= 1
            ecogCAR.data(good_elecs,:) = ecogDS.data(good_elecs,:);
            good_elecs = []; %if there is only one good elec in bank, do not re-reference that channel
        end

        %re-reference (if at least 2 channels are ok in bank)
        for e = 1:length(good_elecs)
            ecogCAR.data(good_elecs(e),:) = ecogDS.data(good_elecs(e),:) - CAR;
        end
    end

    save(CAR_path,'ecogCAR')

end
