function ram_time
% function ram_time
%   Time-domain scrolling plot of data.
% 
%   DR 01/2015

% parameters
ddir = 'G:\ditka\'; % directory
dtank = 'ditka_rec_DT1_032715'; % tank
dblock = 1; % block
ch = 11; % channels
bpf = [0.1 500]; % bandpass filter (Hz)

% constants
fs = 24414.06; % sampling rate (Hz)

% load data
cd([ddir dtank '\Block-' num2str(dblock)]);
[blp,alp] = butter(2,bpf/(fs/2)); 
for ich = 1:length(ch)
    fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(ch(ich)) '.sev'],'r');
    fseek(fid,40,'bof');
    dat = fread(fid,inf,'single');
    fclose(fid);
    dat = filtfilt(blp,alp,dat); % anti-aliasing filter
    dat = downsample(dat,10); % downsample
    if ich == 1, DAT = zeros(length(dat),length(ch)); end
    DAT(:,ich) = dat;
end
fs = round(fs/10);
clear dat;

% interactive plot
figure('Name','oblio','NumberTitle','off','Units','normalized','Position',[1/8 1/8 3/4 3/4]);
ha = axes('Box','off','Units','normalized'); xlabel('sec'); ylabel('\muV'); hold on;
pos = get(ha,'Position');
hs = uicontrol('Style','slider','Units','normalized','Position',[pos(1), pos(2)+pos(4), pos(3), 0.025]);
hc = uicontrol('Style','popupmenu','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.1, 0.05, 0.1],'String',{'20','10','5','2','1','0.5'},'Value',2,'Callback',@PlotCallback);
uicontrol('Style','text','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4), 0.05, 0.02],'String','sec','BackgroundColor',get(gcf,'Color'));
he = uicontrol('Style','edit','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.2, 0.05, 0.025],'String','500','Callback',@PlotCallback);
uicontrol('Style','text','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.175, 0.05, 0.02],'String','uV','BackgroundColor',get(gcf,'Color'));
hb = uicontrol('Style','checkbox','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.25, 0.05, 0.025],'String','CAR','Value',0,'BackgroundColor',get(gcf,'Color'),'Callback',@PlotCallback);
% uicontrol('Style','text','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.175, 0.05, 0.02],'String','uV','BackgroundColor',get(gcf,'Color'));
N = size(DAT,1); % number of samples (minus header)
T = N/fs; % sec
set(hs,'Min',1,'Max',N,'SliderStep',[0.5/T, 2/T],'Value',1,'Callback',@PlotCallback);
gdat = guidata(gcf);
gdat.dat = DAT;
gdat.fs = fs; gdat.T = T;
gdat.hs = hs; gdat.hc = hc;
gdat.he = he; gdat.hb = hb;
guidata(gcf,gdat);
PlotCallback(hs);
end

function PlotCallback(hObject,~)
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
cdat = dat(val:val+ds-1,:);
scale = str2double(get(he,'String'));
car = get(hb,'Value');
if car % common average reference
    cdat = cdat - mean(cdat,2)*ones(1,size(cdat,2));
end

% plot
cla;
t = ((0:size(cdat,1)-1)+val)/fs;
for ich = 1:size(dat,2)
    chdat = cdat(:,ich);
    chdat = (chdat-median(chdat))/scale;
    plot(t,chdat+ich,'k');
end
set(gca,'XLim',[val, val+ds-1]/fs,'YLim',[0 size(dat,2)+1],'YTick',1:size(dat,2));
set(hs,'SliderStep',[0.5*dx/T, 2*dx/T]);
end
