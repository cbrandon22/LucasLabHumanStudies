function ram_time_cage
% function ram_time_cage
%   Time-domain scrolling plot of data collected wirelessly in the cage.
% 
%   DR 03/2015

% parameters
ddir = '/Users/tnl/Desktop/'; % directory
subj = 'ditka'; % subject
tank = 'ditka_cage_051315'; % tank 
fs = 1000; % sampling rate (Hz)
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
DAT = reshape(DAT,nch,[])';

% interactive plot
figure('Name',tank,'NumberTitle','off','Units','normalized','Position',[1/8 1/8 3/4 3/4]);
ha = axes('Box','off','Units','normalized'); xlabel('hr'); ylabel('\muV'); hold on;
pos = get(ha,'Position');
hs = uicontrol('Style','slider','Units','normalized','Position',[pos(1), pos(2)+pos(4), pos(3), 0.025]);
hc = uicontrol('Style','popupmenu','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.1, 0.05, 0.1],'String',{'20','10','5','3','2','1','0.5'},'Value',3,'Callback',@PlotCallback);
uicontrol('Style','text','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4), 0.05, 0.02],'String','sec','BackgroundColor',get(gcf,'Color'));
he = uicontrol('Style','edit','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.2, 0.05, 0.025],'String','500','Callback',@PlotCallback);
uicontrol('Style','text','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.175, 0.05, 0.02],'String','uV','BackgroundColor',get(gcf,'Color'));
hb = uicontrol('Style','checkbox','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.25, 0.05, 0.025],'String','CAR','Value',0,'BackgroundColor',get(gcf,'Color'),'Callback',@PlotCallback);
hp = uicontrol('Style','pushbutton','Units','normalized','Position',[pos(1)+pos(3), pos(2), 0.07, 0.1],'String','PSD','Callback',@PlotCallback,'Enable','off');
N = size(DAT,1); % number of samples
T = N/fs; % sec
set(hs,'Min',1,'Max',N,'SliderStep',[0.5/T, 2/T],'Value',1,'Callback',@PlotCallback);
gdat = guidata(gcf);
gdat.dat = DAT;
gdat.fs = fs; gdat.T = T;
gdat.hs = hs; gdat.hc = hc;
gdat.he = he; gdat.hb = hb; gdat.hp = hp;
guidata(gcf,gdat);
PlotCallback(hs);
end

function PlotCallback(hObject,~)
% load data
gdat = guidata(gcf);
dat = gdat.dat;
fs = gdat.fs; T = gdat.T;
hs = gdat.hs; hc = gdat.hc;
he = gdat.he; hb = gdat.hb; hp = gdat.hp;
val = round(get(hs,'Value')); % start sample
dxs = get(hc,'String');
dx = str2double(dxs{get(hc,'Value')});
if dx > 2, set(hp,'Enable','on'); else set(hp,'Enable','off'); end
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
    plot(t/3600,chdat+ich,'k');
end
set(gca,'XLim',[val, val+ds-1]/fs/3600,'YLim',[0 size(dat,2)+1],'YTick',1:size(dat,2));
set(hs,'SliderStep',[0.5*dx/T, 2*dx/T]);
if hObject==hp
    figure('Name',[num2str(val/fs,'%4.2f') '-' num2str((val+ds-1)/fs,'%4.2f')],'NumberTitle','off','Units','normalized','Position',[1/2 1/3 1/3 1/3]);
    for ich = 1:size(dat,2)
        chdat = cdat(:,ich)-median(cdat(:,ich));
        [P,f] = pwelch(chdat,2*fix(fs),fix(fs),2*fix(fs),fs); % 2-sec windows, 1-sec overlap
%         [P,f] = pmtmwelch(chdat,3,2*fix(fs),fix(fs),2*fix(fs),fs);
        ind = find(f>=.1 & f<=200);
        plot(f(ind),P(ind),'k'); hold on;
    end
    axis tight; set(gca,'Box','off','XScale','log','YScale','log','YLim',[1e-1 1e5]);
    xlabel('Hz'); ylabel('\muV^2/Hz');
end
end
