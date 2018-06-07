function EMUstim
% function EMUstim
%   Look at neural recordings around stimulation events delivered to
%   humans in EMU.
%
%   DR 03/2018

% parameters
ddir = '/Volumes/HumanStudies/HumanStudies/corticalStim/eeg/HUP128/processed'; % data directory
sind = 11; % stimulus index
win = [-15 30]; % peristimulus time window (sec)
scale = 2000; % scale factor (uV?)

% constants
fs = 1024; % sampling rate (Hz)

% data
cd(ddir);
[tstim,lead,elec1,elec2] = stim_info(sind);
load([lead '_dat.mat'],'-mat');
t = (1:size(eeg,2))/fs;
[~,istim] = min(abs(t-tstim));
swin = (round(win(1)*fs):round(win(2)*fs))+istim;
datwin = eeg(:,swin);
datwin(elec1,:) = NaN*ones(1,size(datwin,2));
datwin(elec2,:) = NaN*ones(1,size(datwin,2));
twin = (swin-istim)/fs;

% plot
figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/8 1/8 3/4 3/4],'Color','w');
cmap = colormap(lines(size(eeg,1)));
for ii = 1:size(eeg,1)
    plot(twin,datwin(ii,:)/scale+ii,'Color',cmap(ii,:),'Tag',num2str(ii)); hold on;
end
axis tight; set(gca,'Box','off','YLim',[0 size(datwin,1)+1],'XGrid','on');
xlabel('time (s)'); ylabel('electrode')
end

function [tstim,lead,elec1,elec2] = stim_info(ind)
%   Returns info for specified stim event

start_time = '15:28:14';%  		Beginning of Recording

stim_time = {'16:01:00','LA',6,7;%1  		Closed relay to LA6 and LA7
            '16:03:42','LA',6,7;%2  		Closed relay to LA6 and LA7
            '16:04:17','LA',6,7;%3  		Closed relay to LA6 and LA7
            '16:04:56','LA',6,7;%4  		Closed relay to LA6 and LA7
            '16:06:24','LB',5,6;%5  		Closed relay to LB5 and LB6
            '16:10:50','LB',7,8;%6  		Closed relay to LB7 and LB8
            '16:13:21','LC',5,6;%7  		Closed relay to LC5 and LC6
            '16:14:54','LC',5,6;%8  		Closed relay to LC5 and LC6
            '16:15:50','LC',5,6;%9  		Closed relay to LC5 and LC6
            '16:22:36','LC',7,8;%10  		Closed relay to LC7 and LC8
            '16:25:19','LX',4,5;%11  		Closed relay to LX4 and LX5
            '16:28:20','LX',5,6;%12  		Closed relay to LX5 and LX6
            '16:31:08','LW',3,4;%13  		Closed relay to LW3 and LW4
            '16:32:32','LZ',5,6;%14  		Closed relay to LZ5 and LZ6
            '16:34:33','LZ',7,8;%15  		Closed relay to LZ7 and LZ8
            '16:35:41','LY',5,6;%16  		Closed relay to LY5 and LY6
            '16:40:23','LY',7,8;%17  		Closed relay to LY7 and LY8
            '16:53:56','LY',6,7;%18  		Closed relay to LY6 and LY7 unable to name words with stim
            '17:08:35','LW',3,4;%19  		Closed relay to LW3 and LW4
            '17:09:28','LZ',6,7;%20  		Closed relay to LZ6 and LZ7
            '17:12:07','LX',5,6;%21  		Closed relay to LX5 and LX6
            '17:13:18','LY',6,7};%22  		Closed relay to LY6 and LY7 disrupt speech comprehension

t1 = datevec(start_time,'HH:MM:SS');
t2 = datevec(stim_time{ind,1},'HH:MM:SS');
tstim = etime(t2,t1);
lead = stim_time{ind,2};
elec1 = stim_time{ind,3};
elec2 = stim_time{ind,4};
end