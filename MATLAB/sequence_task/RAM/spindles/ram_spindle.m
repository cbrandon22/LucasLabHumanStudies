function [ep,t] = ram_spindle(varargin)
% function [ep,t] = ram_pepisode(varargin)
%   Oscillatory episode detection using Pepisode method (Caplan et al
%   2001).
%
%   DR 03/2015

% parameters
ddir = '/Users/tnl/Desktop/'; % directory
subj = 'ditka'; % subject
tank = 'ditka_cage_040615'; % tank
% block = 1; % block
nch = 1; % number of channels
dns = 50; % downsample factor
wavl = 'morl'; % wavelet
fev = 100:2:200; % evaluation frequencies
tcyc = 3; % duration threshold (cycles)
fs = 1000; % Hz
winsz = 5; % window size for spectral estimates (s) - determines frequency resolution
if nargin
    v2struct(varargin{1});
end

% data
cd([ddir subj filesep tank]);
    disp('loading data...');
    Nfl = length(dir('*_*.txt'));
    dat = [];
    for ii = 0:Nfl-1
        fnm = [tank(end-5:end-2) '_' num2str(ii) '.txt'];
        fid = fopen(fnm,'r');
        dat = [dat; fread(fid,inf,'uint16')];
        fclose(fid);
    end
    lastfullcycle = floor(length(dat)/nch)*nch;
    dat(lastfullcycle+1:end) = [];
    dat = reshape(dat,nch,[]);
    N = size(dat,1); % samples
    S = fix(N/fs)-(winsz-1); % seconds


if dns>1
    [b,a] = butter(2,round(fs/(2.5*dns))/(fs/2)); % antialiasing filter
    dat = filtfilt(b,a,dat);
end
dat = downsample(dat,dns);
fs = round(fs/dns);
t = (1:length(dat))/fs;
save([ddir subj filesep 'processed' filesep tank '_vt1'],'dat','t','N','S','fs');