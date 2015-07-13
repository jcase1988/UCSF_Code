%% Retrieve subject globals (e.g., dtdir, anadir, figdir, anadir_day, fig_dir_day)
%  from subj_globals file in parent directory

function get_subj_globals(subj, block)
    
    % Load the single subj_globals file with all subj/block data
    load([UCSF_dir filesep 'subj_globals.mat']);
    
    % Retrieve all variable names associated with subj/block (dtdir,
    % anadir etc.)
    fields = fieldnames(subj_globals.(subj).(block));
    
    % Return all variables to caller workspace
    for i = 1:length(fields)
        assignin('caller',fields{i},subj_globals.(subj).(block).(fields{i}));
    end
        
end