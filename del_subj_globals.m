function del_subj_globals(subj,block,varargin)

    %UCSF_dir defined in function UCSF_dir()
    parentdir = UCSF_dir();
    sep = filesep;

    %Load subject globals
    load([parentdir sep 'subj_globals.mat']);
    
    for iarg = 1:length(varargin)        
        subj_globals.(subj).(block) = rmfield(subj_globals.(subj).(block),varargin{iarg});
    end
    
    %Save subject globals
    sav = input('Delete these data? Yes (1) or No (2) ');
    if sav == 1
        save([parentdir sep 'subj_globals.mat'],'subj_globals');
    end
    
end