function create_subj_globals(subj)
    parentdir = UCSF_dir;
    sep = filesep;
    subjdir = [UCSF_dir sep subj sep];
    
    %get contents of subj directory
    d = dir(subjdir);
    d = regexpi({d.name},'.*B\d.*','match');
    block_names = [d{:}];   
   
    
    dtdir = ['/Users/johncase/Documents/UCSF Data/' subj '/' block '/'];
    load([dtdir subj '_' block])
    srate = ecogDS.sampFreq;
end