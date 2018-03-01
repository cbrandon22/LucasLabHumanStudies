function ditka_pac_pop
% function ditka_pac_pop
%
%   DR 03/2015

% parameters
p.ddir = 'D:\'; % directory
p.subj = 'ditka'; % subject
tanks = [...%'ditka_rec_DT1_031115';...% data tanks
         ...%'ditka_rec_DT1_031715';...
         ...'ditka_rec_DT1_031915';...
         'ditka_rec_DT1_032315'];         
p.block = 1; % data block
p.ch = 5; % channel(s)
p.fP = [3 7]; % bandpass cutoffs - phase range
fAs = [35 45; 60 90]; % bandpass cutoffs - amplitude range
p.delta = pi/8; % phase bin size (rad)
p.cph = 1; % correct for waveform asymmetry? (0 or 1)

% preferred phase
Nt = size(tanks,1);
Na = size(fAs,1);
PP = zeros(Nt,Na);
for ii = 1:Nt
    p.tank = tanks(ii,:);
    for jj = 1:Na
        p.fA = fAs(jj,:);
        [~,~,pp] = ram_pac(p);
        PP(ii,jj) = pp(2); % preferred phase at high-theta amplitudes
    end
end
% close all;

% plot
figure('Name','ditka pac pop','NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
x = 0:pi/64:2*pi;
plot(x*180/pi,sin(x),'k','LineWidth',2); hold on;
axis tight; set(gca,'Box','off','Visible','off');
clr = ['b' 'r'];
PP = PP+pi/2; % shift phases from cosine to sine
for ii = 1:Nt
    for jj = 1:Na
        h(jj) = plot(PP(ii,jj)*180/pi,sin(PP(ii,jj)),'Color',clr(jj),'Marker','o','MarkerSize',14,'MarkerFaceColor',clr(jj),'Tag',[num2str(ii) '-' num2str(jj)]);
    end
end
legend(h,'low-gamma','high-gamma');
legend('boxoff');
suptitle('preferred theta phase');