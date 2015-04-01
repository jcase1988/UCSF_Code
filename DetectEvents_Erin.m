
function [eventtimes,eventnames,aga,ushu,ifi,clicks]=DetectEvents_Erin(dtpath,wpath,anin_to_use,names,tasktype)
%function [eventtimes,eventnames,trialstims,faces,phrases,clicks]=DetectEvents_Erin(dtpath,wpath,anin_to_use,names,tasktype);
%function [eventtimes,eventnames]=DetectEvents_Erin(dtpath,wpath,anin_to_use,names)


%dtpath = ANIN data path  dtpath='/Users/Erin/Documents/MATLAB/' %path to
%ANIN.htk files
%wpath = path to sample wav files wpath='/Users/Erin/Documents/MATLAB/ecog
%preprocessing/Erin/EmotIDsounds/block6' %path to wave files
%anin_to_use = ANIN channel %use 2 or 3 if possible
%names = cell array of names of sample wav files

new_fs = 20000;
corr_thresh=0.3;
evind = 1;
all_conf=[];
eventnames=cell(1,1);
eventtimes=[];
sentence=cell(length(names),1);

%keyboard
wname = sprintf('%sANIN%d.htk', dtpath, anin_to_use);
%keyboard
[anin2,anin2_fs] = readhtk(wname); %load analog channel

anin2 = resample(anin2,round(new_fs),round(anin2_fs)); %resample analog signals to 16kHz
anin2 = anin2-mean(anin2); % subtract mean to get mean 0
anin2 = 0.99*anin2/max(abs(anin2));
anin2=-anin2;

for stimnum=1:length(names) 
    wname = [wpath filesep names{stimnum}];
    [sentence{stimnum}, sent_fs] = wavread(wname);
    sentence{stimnum} = resample(sentence{stimnum}, new_fs, round(sent_fs)); 
    sentence{stimnum}=sentence{stimnum}(:,1);
end 
%keyboard
for stimnum=1:length(names)

    % Find where this sentence audio occurs in ANIN2
    fftlen=length(anin2);
    match_filter= cconv(flipud(sentence{stimnum}), anin2(:), fftlen); % THIS IS MUCH FASTER THAN CONV!
  
    
    c=1;
    while c>corr_thresh  
        % sort by maximum of the convolution, this tells you where the END
        % of the events occurred. We sort in case there is more than one
        % example of the sentence detected in the TDT signal
        [~, maxind]=max(match_filter);
        cc=[];
        xc=[];
        for altstims=1:length(names)
            [stimnum altstims]
            start_time=maxind-length(sentence{altstims})+1; % this is where the events start
            end_time=(start_time+length(sentence{altstims})-1);    
            matched_sentence_segment = anin2(start_time:end_time);
            %if stimnum==altstims
            %    xc=corrcoef(sentence{altstims}, matched_sentence_segment);
            %    xc=xc(1,2);
            %else
                xc=xcorr(sentence{altstims}, matched_sentence_segment,floor(length(sentence{altstims})/2),'coeff'); % correlation between sentence and the "match"
            %end
            cc(altstims)=max(xc);
            %all_conf = [all_conf; cc(1,2)];
        end
        if cc(stimnum)==max(cc)
            start_time=maxind-length(sentence{stimnum})+1;
            end_time=(start_time+length(sentence{stimnum})-1); 
            eventtimes(evind,1)=start_time/new_fs;
            eventtimes(evind,2)=end_time/new_fs;
            eventnames{evind,1}=names{stimnum};
            evind=evind+1;
            
        end
        match_filter(floor(start_time+length(sentence{stimnum})/2):ceil(end_time+length(sentence{stimnum})/2)) = 0;
        anin2(start_time:end_time) = 0;
        c=cc(stimnum);
    end
end

[eventtimes evind]=sortrows(eventtimes);
eventnames=eventnames(evind);


if tasktype==1
clickcyc=0;
facecyc=0;
soundstimcyc=0;
trialstimcyc=0;
for k=1:length(eventnames)
    if strcmpi(eventnames{k},'click')==1
                clickcyc=clickcyc+1;
        clicks(clickcyc,:)=eventtimes(k,:);
    elseif strcmpi(eventnames{k},'zoop')==1
        facecyc=facecyc+1;
        trialstimcyc=trialstimcyc+1;
        faces(facecyc,:)=eventtimes(k,:);
        trialstims(trialstimcyc,:)=eventtimes(k,:);
    else
        soundstimcyc=soundstimcyc+1;
        trialstimcyc=trialstimcyc+1;
        phrases(soundstimcyc,:)=eventtimes(k,:);
        trialstims(trialstimcyc,:)=eventtimes(k,:);
    end
end

elseif tasktype==2
agacyc=0;
ificyc=0;
ushucyc=0;
ipicyc=0;
trialstimcyc=0;
clickcyc=0;
for k=1:length(eventnames)
    if strcmpi(eventnames{k},'aagaa')==1
        agacyc=agacyc+1;
        aga(agacyc,:)=eventtimes(k,:);
        trialstimcyc=trialstimcyc+1;
        trialstims(trialstimcyc,:)=eventtimes(k,:);
    elseif strcmpi(eventnames{k},'iyfiy')==1
        ificyc=ificyc+1;
        trialstimcyc=trialstimcyc+1;
        ifi(ificyc,:)=eventtimes(k,:);
        trialstims(trialstimcyc,:)=eventtimes(k,:);
    elseif strcmpi(eventnames{k},'ihpih')==1
        ipicyc=ipicyc+1;
        trialstimcyc=trialstimcyc+1;
        ipi(ipicyc,:)=eventtimes(k,:);
        trialstims(trialstimcyc,:)=eventtimes(k,:);
    elseif strcmpi(eventnames{k},'click')==1
        clickcyc=clickcyc+1;
        clicks(clickcyc,:)=eventtimes(k,:);
    else
        ushucyc=ushucyc+1;
        trialstimcyc=trialstimcyc+1;
        ushu(ushucyc,:)=eventtimes(k,:);
        trialstims(trialstimcyc,:)=eventtimes(k,:);
    end
end
end
        
        