%% Initialize subj_globals: add_subj_globals('EC77','B26')
%% or add more fields:      add_subj_globals('EC77','B26','Day1',bad_epochs)

function add_subj_globals(subj,block,varargin)
        
    %UCSF_dir defined in function UCSF_dir()    
    parentdir = UCSF_dir();
    sep = filesep;
    
    %Load subject globals
    load([parentdir sep 'subj_globals.mat']);  
    
    %If you are not assigning anything beyond file paths
    if isempty(varargin)
        
        %If 'All' 
        if strcmp('All',block)
            figdir_day = [parentdir sep subj sep 'All' sep 'figures' sep];    
            anadir_day = [parentdir sep subj sep 'All' sep 'analysis' sep];
            subj_globals.(subj).All.figdir_all = figdir_day;
            subj_globals.(subj).All.anadir_all = anadir_day;
        else

            %Determine file paths
            subjdir = [UCSF_dir sep subj sep];    
            dtdir = [parentdir sep subj sep subj block sep 'data' sep];
            figdir = [parentdir sep subj sep subj block sep 'figures' sep];        
            anadir = [parentdir sep subj sep subj block sep 'analysis' sep];

            load([dtdir subj '_' block '.mat'])
            srate = ecogDS.sampFreq;

            %Assign to subj_globals.mat

            subj_globals.(subj).(block).srate = srate;
            subj_globals.(subj).(block).subjdir = subjdir;
            subj_globals.(subj).(block).dtdir = dtdir;
            subj_globals.(subj).(block).anadir = anadir;   
            subj_globals.(subj).(block).figdir = figdir;  
        end

    %if assigning field names to subj_globals
    else 
        
        
        for iarg = 1:length(varargin)
            
            %check to see if argment is "Day"
            low_arg = lower(inputname(2+iarg));
            if strcmp('day',low_arg(1:3))
                figdir_day = [parentdir sep subj sep varargin{iarg} sep 'figures' sep];    
                anadir_day = [parentdir sep subj sep varargin{iarg} sep 'analysis' sep];
                subj_globals.(subj).(block).figdir_day = figdir_day;
                subj_globals.(subj).(block).anadir_day = anadir_day;
                subj_globals.(subj).(block).day = varargin{iarg};
            
                %Create Day subj globals 
                subj_globals.(subj).(varargin{iarg}).figdir_day = figdir_day;
                subj_globals.(subj).(varargin{iarg}).anadir_day = anadir_day;
                
                if isfield(subj_globals.(subj).(varargin{iarg}),'blocks')
                    subj_globals.(subj).(varargin{iarg}).blocks = sort(unique([subj_globals.(subj).(varargin{iarg}).blocks block]));
                else
                    subj_globals.(subj).(varargin{iarg}).blocks = {block};
                end
            % 'All' blocks
            elseif strcmp('all',low_arg(1:3)) && ~strcmp('allstimtimes',low_arg)
                figdir_day = [parentdir sep subj sep 'All' sep 'figures' sep];    
                anadir_day = [parentdir sep subj sep 'All' sep 'analysis' sep];
                subj_globals.(subj).All.figdir_all = figdir_day;
                subj_globals.(subj).All.anadir_all = anadir_day;
            
            % if not a "Day" or 'All' variable
            else
                vstring = inputname(2+iarg); 
                subj_globals.(subj).(block).(vstring) = varargin{iarg};
            end
        end
    
    end
    
    save([parentdir sep 'subj_globals.mat'],'subj_globals');  
    
end