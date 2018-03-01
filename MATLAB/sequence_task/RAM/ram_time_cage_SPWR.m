function ram_time_cage_SPWR
% function ram_time_cage_SPWR
%   Time-domain scrolling plot of data collected wirelessly in the cage:
%   version to enable visualization of sharp wave-ripples.
% 
%   DR 08/2016

% parameters
ddir = '/Users/tnl/Desktop/'; % directory
subj = 'ditka'; % subject
tank = 'ditka_cage_040615'; % tank 
fs = 1000; % sampling rate (Hz)
nch = 1; % number of channels
ch = 1; % channel to analyze

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
DAT = reshape(DAT,nch,[])';
DAT = DAT(:,ch);
[bB,aB] = butter(2,[150 250]/(fs/2)); % ripple bandpass filter
RIP = filtfilt(bB,aB,DAT); % ripple-band activity

% interactive plot
figure('Name',tank,'NumberTitle','off','Units','normalized','Position',[1/8 1/8 3/4 3/4]);
ha = axes('Box','off','Units','normalized'); xlabel('hr'); ylabel('\muV'); hold on;
pos = get(ha,'Position');
hs = uicontrol('Style','slider','Units','normalized','Position',[pos(1), pos(2)+pos(4), pos(3), 0.025]);
hc = uicontrol('Style','popupmenu','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.1, 0.05, 0.1],'String',{'20','10','5','3','2','1','0.5'},'Value',3,'Callback',@PlotCallback);
uicontrol('Style','text','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4), 0.05, 0.02],'String','sec','BackgroundColor',get(gcf,'Color'));
he = uicontrol('Style','edit','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.2, 0.05, 0.025],'String','1000','Callback',@PlotCallback);
uicontrol('Style','text','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.175, 0.05, 0.02],'String','wideband (uV)','BackgroundColor',get(gcf,'Color'));
hf = uicontrol('Style','edit','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.3, 0.05, 0.025],'String','50','Callback',@PlotCallback);
uicontrol('Style','text','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.275, 0.05, 0.02],'String','ripple (uV)','BackgroundColor',get(gcf,'Color'));
N = size(DAT,1); % number of samples
T = N/fs; % sec
set(hs,'Min',1,'Max',N,'SliderStep',[0.5/T, 2/T],'Value',1,'Callback',@PlotCallback);
gdat = guidata(gcf);
gdat.dat = DAT; gdat.rip = RIP;
gdat.fs = fs; gdat.T = T;
gdat.hs = hs; gdat.hc = hc;
gdat.he = he; gdat.hf = hf;
guidata(gcf,gdat);
PlotCallback(hs);
end

function PlotCallback(~,~)
% load data
gdat = guidata(gcf);
dat = gdat.dat; rip = gdat.rip;
fs = gdat.fs; T = gdat.T;
hs = gdat.hs; hc = gdat.hc;
he = gdat.he; hf = gdat.hf;
val = round(get(hs,'Value')); % start sample
dxs = get(hc,'String');
dx = str2double(dxs{get(hc,'Value')});
ds = fix(fs*dx);
cdat = dat(val:val+ds-1);
crip = rip(val:val+ds-1);
scaleE = str2double(get(he,'String'));
scaleF = str2double(get(hf,'String'));

% plot
cla;
t = ((0:size(cdat,1)-1)+val)/fs;
plot(t/3600,(cdat-median(cdat))/scaleE+2,'k');
plot(t/3600,(crip-median(crip))/scaleF+1,'k');
set(gca,'XLim',[val, val+ds-1]/fs/3600,'YLim',[0 3],'YTick',1:2,'YTickLabel',{'ripple','lfp'});
set(hs,'SliderStep',[0.5*dx/T, 2*dx/T]);
end
