function ram_saccade(varargin)
% function ram_saccade(varargin)
%   Find indices of saccades.
%
%   DR 08/2015

% parameters
ddir = 'D:\'; % data directory
subj = 'fidel'; % subject
dtank = 'fidel_VPLT_DT1_111715'; % data tank
dblock = 1; % data block
if nargin
    v2struct(varargin{1});
end

% constants
fs = 24414.06; % sampling rate (Hz)

% load eye tracker channels
cd([ddir subj]);
cd([dtank '\Block-' num2str(dblock)]);
fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(17) '.sev'],'r');
fseek(fid,40,'bof');
eye0 = fread(fid,inf,'single')/1e6; % eye X
fclose(fid);
fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(18) '.sev'],'r');
fseek(fid,40,'bof');
eye1 = fread(fid,inf,'single')/1e6; % eye Y
fclose(fid);

% find saccades
[Sindex, Ssize, Sdir] = saccadeFinder(eye0,eye1,fs);
saccade = [Sindex Ssize Sdir]; % N x [onset offset size direction]

% load processed data
try load([ddir subj '\processed\' dtank '_' num2str(dblock) '.mat'],'-mat');
catch disp('run ''ram_stimtrig'' first'); return; end

% save
save([ddir subj '\processed\' dtank '_' num2str(dblock) '.mat'],'saccade','-append','-mat');