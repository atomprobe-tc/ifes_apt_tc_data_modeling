raw = fopen(["/home/kaiobach/Research/paper_paper_paper/joss_nomad_apt/bb_analysis/data/fra_rouen_karam/Experimental Analysis of LaB6 ; Negative Pulse ; DC Voltage 6.6 kV ; Amplitudes 2.5 kV and 3 kV.ato"]);
attoread = inf; %max number of atoms
data.version = fread(raw,2,'uint32','l');
display(data.version(2));
header = 5000;
oneat = 320;
fseek(raw,header,'bof');
wpx = fread(raw,attoread,'int16',oneat/8-2,'b');
fseek(raw,header+2,'bof');
wpy = fread(raw,attoread,'int16',oneat/8-2,'b');
fclose(raw);
