clear

subj = 'EC71';
block = 'B11';

fp = 1:0.1:15.1;
fp_bandwidth = 0.5:0.1:5.1;


get_subj_globals(subj,block);
grid_loc = [64:-8:1 63:-8:1 62:-8:1 61:-8:1 60:-8:1 59:-8:1 58:-8:1 57:-8:1];

phase_elec = 1:48;
amp_elec = 90;

fgrid = figure; 
set(fgrid, 'Position', [1 1 2400 1260]);
for iElec = phase_elec
    load([anadir 'MI/e' num2str(iElec) '_e' num2str(amp_elec) '.mat'])
    
    pos = find(grid_loc == iElec);
    subplot(8,8,pos)
    heatmap(MI,fp,fp_bandwidth);
end
    