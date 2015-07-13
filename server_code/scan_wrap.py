import sys
from func_scan_frequencies_for_CFC import scan_freq

subj = 'EC82'
block = 'B93'

phase_elec = int(sys.argv[1])
amp_elec = int(sys.argv[2])

print('Phase_elec ' + str(phase_elec) + ', Amp_elec ' + str(amp_elec)) 
scan_freq(subj,block,phase_elec,amp_elec)
