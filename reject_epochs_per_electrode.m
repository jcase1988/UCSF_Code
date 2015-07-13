function reject_epochs_per_electrode(subj,block,depth)

%depth: 0 = all electrodes
%       1 = depth electrodes only

get_subj_globals(subj,block)
load([dtdir subj '_' block '_CAR.mat'])

if ~exist('per_chan_bad_epochs')
    per_chan_bad_epochs = cell(size(ecogCAR.data,1),1);
end

if depth
    if strcmp(subj,'EC71')
        ecogCAR.data(1:128,:) = 0;
    else
        ecogCAR.data(1:256,:) = 0;
    end
end

%Define event structure
events = allstimtimes;

for i = 1:length(events)
    event_struct(i).type='aud';
    event_struct(i).latency=round(events(i,1)*srate);
    event_struct(i).duration=round( (events(i,2)-events(i,1)) *srate);
end

%Filter data into high gamma envelop 
ecogCAR = hgf(ecogCAR);
eegplot(ecogCAR.data,'srate',ecogCAR.sampFreq,'events',event_struct)


%User input data
inp = '';
all_inps = {};
while true
    inp = input('','s');
    
    if strcmp(inp,'q')
        break;
    end
    
    all_inps = [all_inps inp];
    
end

%Parse input data
for iInp = 1:length(all_inps)
    inp = all_inps{iInp};
    
    %if set of electrodes (e.g., "[2 4 6] 100 101")
    if ~isempty(strfind(inp,'[')) 
        beg = strfind(inp,'[')+1;
        en = strfind(inp,']')-1;
        elec_split = strsplit(inp(beg:en),' ');
        epoch_split = strsplit(strtrim(inp(en+2:end)), ' ');
        beg_epoch = str2num(epoch_split{1});
        end_epoch = str2num(epoch_split{2});
                
        for iElec = 1:length(elec_split)
            
            %if there is a ":" in a bracket
            if ~isempty(strfind(elec_split{iElec},':'))
                elec_range = strsplit(elec_split{iElec},':');
                beg_elec = str2num(elec_range{1});
                end_elec = str2num(elec_range{2});
                
                for jElec = beg_elec:end_elec
                    per_chan_bad_epochs{jElec} = [per_chan_bad_epochs{jElec} ; beg_epoch end_epoch];
                end
            else
                elec = str2num(elec_split{iElec});
                per_chan_bad_epochs{elec} = [per_chan_bad_epochs{elec} ; beg_epoch end_epoch];
            end
        end
    %if just one electrode or a range without brackets 
    else
        inp_split = strsplit(inp,' ');    
        elecs = str2num(inp_split{1});
        beg_epoch = str2num(inp_split{2});
        end_epoch = str2num(inp_split{3});

        for elec = elecs
            per_chan_bad_epochs{elec} = [per_chan_bad_epochs{elec} ; beg_epoch end_epoch];
        end
    end
end

%Save bad epochs
sav = input('Save bad epoch? Yes (1) or No (2) ');
if sav == 1
    add_subj_globals(subj,block,per_chan_bad_epochs);
end
    