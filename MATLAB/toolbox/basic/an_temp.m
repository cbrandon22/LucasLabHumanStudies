%temp
bpPairs = [1 2;1 9;2 3;2 10;3 4;3 11;4 5;4 12;5 6;5 13;6 7;6 14;7 8;7 15;9 10;10 11;11 12;12 13;13 14;14 15;15 16;8 16];
for i = 1:length(bpPairs)
    hold all
    figure;
    samplingRate = 
    endTime = (events(length(events)).lfpoffset)
    lfp = an_getlfp_ms_wrapper(bpPairs(i,:),events(5),1200000,-100000,0);
    %plotSEM(1:1200000,nanmean(lfp,1),SEM(lfp))
    xlabel('Time from Sound (ms)')
    ylabel ('voltage (mv)')
    set(gca,'ylim',[-100 100])
    plot([0 0],get(gca,'ylim'))
    plot(lfp)
    pause
    
end

function [dat,samplerate] = loadData(filename)
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