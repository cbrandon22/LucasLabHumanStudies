function [dat,samplerate] = loadeData(filename)
%LOADDATA - Load EEG data from a file
%

% get the data format
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(filename);

% Open and load the file
fid = fopen(filename, 'r','l');
dat =  fread(fid, inf, dataformat);
fclose(fid);

% apply the gain
dat = dat.*gain;
end