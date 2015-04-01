%% Return the top-most directory where the UCSF data is stored.
% This is the only function that should need to be modified across
% platforms (hopefully)
% Used for create_subj_globals.m
function [d] = UCSF_dir()
    d = '~/Documents/UCSF Data';
    return