function nlx_time
% function nlx_time
% 	Scrolling time-domain plot of one raw signal recorded on Neuralynx
% 	system.
%
%   DR 10/2017

% parameters
dirs = le_dirs('oddball');
ddir = fullfile(dirs.data,'eeg/');
%ddir = 'E:\TNL_Data\oddball\eeg\';  % data directory
subj = 'HUP149_e';                  % subject
elec = 16;                          % electrode

% load data
cd([ddir subj '/raw']);
fl = dir('*-*');
cd(fl(1).name);
[data,info,~,~] = load_ncs2(['CSC' num2str(elec) '.ncs']);
data = double(data);
fs = info.actualSampleRate;

% interactive plot
figure('Name',subj,'NumberTitle','off','Units','normalized','Position',[1/8 1/8 3/4 3/4]);
ha = axes('Box','off','Units','normalized'); xlabel('sec'); ylabel('\muV'); hold on;
pos = get(ha,'Position');
hs = uicontrol('Style','slider','Units','normalized','Position',[pos(1), pos(2)+pos(4), pos(3), 0.025]);
hc = uicontrol('Style','popupmenu','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.1, 0.05, 0.1],'String',{'20','10','5','2','1','0.5'},'Value',2,'Callback',@PlotCallback);
uicontrol('Style','text','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4), 0.05, 0.02],'String','sec','BackgroundColor',get(gcf,'Color'));
he = uicontrol('Style','edit','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.2, 0.05, 0.025],'String','500','Callback',@PlotCallback);
uicontrol('Style','text','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.175, 0.05, 0.02],'String','uV','BackgroundColor',get(gcf,'Color'));
hb = uicontrol('Style','checkbox','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.25, 0.05, 0.025],'String','CAR','Value',0,'BackgroundColor',get(gcf,'Color'),'Callback',@PlotCallback);
N = length(data); % number of samples (minus header)
T = N/fs; % sec
set(hs,'Min',1,'Max',N,'SliderStep',[0.5/T, 2/T],'Value',1,'Callback',@PlotCallback);
gdat = guidata(gcf);
gdat.dat = data;
gdat.fs = fs; gdat.T = T;
gdat.hs = hs; gdat.hc = hc;
gdat.he = he; gdat.hb = hb;
guidata(gcf,gdat);
PlotCallback(hs);
end

function PlotCallback(~,~)
% load data
gdat = guidata(gcf);
dat = gdat.dat;
fs = gdat.fs; T = gdat.T;
hs = gdat.hs; hc = gdat.hc;
he = gdat.he; hb = gdat.hb;
val = round(get(hs,'Value')); % start sample
dxs = get(hc,'String');
dx = str2double(dxs{get(hc,'Value')});
ds = fix(fs*dx);
cdat = dat(val:val+ds-1);

% plot
cla;
t = ((0:length(cdat)-1)+val)/fs;
plot(t,cdat,'k');
set(gca,'XLim',[val, val+ds-1]/fs);
set(hs,'SliderStep',[0.5*dx/T, 2*dx/T]);
end
