function [DAT,t]= ExtractRaw;

    ddir = '/Users/tnl/Desktop/'; % directory
    subj = 'ditka'; % subject
    tank = 'ditka_cage_050515'; % tank 
%     fs = 1000; % sampling rate (Hz)
    nch = 1; % number of channels

    % load data
    cd([ddir subj filesep tank]);
    Nfl = length(dir('*_*.txt'));
    DAT = [];
    for ii = 0:Nfl-1
        fnm = [tank(end-5:end-2) '_' num2str(ii) '.txt'];
        fid = fopen(fnm,'r');
        DAT = [DAT; fread(fid,inf,'uint16')];
        fclose(fid);
    end
    lastfullcycle = floor(length(DAT)/nch)*nch;
    DAT(lastfullcycle+1:end) = [];
    dat = reshape(DAT,nch,[])';
    t = 1:length(dat);
 save([ddir subj filesep 'processed' filesep tank '_vtRaw'],'-v7.3','dat','t');
end