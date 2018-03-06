function ram_stx_cage
% function ram_stx_cage
%   Stacked plot of every stimulus-triggered trial - cage experiment.
%
%   DR 05/2015

% parameters
ddir = 'F:\'; % directory
subj = 'ditka'; % subject
tank = 'ditka_cage_051815'; % tank 
scale = 300; % scale (uV)
xlim = [-10 50]; % (ms)

% data
dscale = 3.3/2^12/1000*1e6; % amplitude scale factor to convert to uV
try 
    load([ddir subj '\processed\' tank '.mat'],'-mat');
    N = size(STA,1);
    STA = dscale*(STA-median(STA(:,t<0),2)*ones(1,length(t))); % center and scale
catch
    error('must run ''ram_sta_cage.m'' first');
end

% plot
figure('Name',tank,'NumberTitle','off','Units','normalized','Position',[1/3 1/8 1/3 3/4],'Color','w');
imagesc([t(1) t(end)],[1 N],STA); hold on;
set(gca,'Box','off','YDir','normal','XLim',xlim,'YLim',[1 N],'TickDir','out');
for ii = 1:N
    plot(t,STA(ii,:)/scale+ii,'Color','k','Tag',num2str(ii)); hold on;
end
 xlabel('time (ms)'); ylabel('stimulus number');
line(xlim(2)*ones(1,2),[0 10],'Color','k');
text(xlim(2),5,[num2str(scale*10) '\muV'],'HorizontalAlignment','left','VerticalAlignment','middle');
