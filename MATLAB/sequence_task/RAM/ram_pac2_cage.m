function ram_pac2_cage
% function ram_pac2_cage
%   Compute phase-amplitude coupling from recording files (no stim):
%   comodulogram plot
% 
%   DR 03/2015

% parameters
ddir = '/Users/tnl/Desktop/'; % directory
subj = 'ditka'; % subject
tank = 'ditka_cage_040615'; % tank 
flnm = '16.08.05_record.txt'; % file
dhr = [8 10]; % time window (hours)
delta = pi/8; % phase bin size (rad)
cph = 1; % correct for waveform asymmetry? (0 or 1)

% constants
dfs = 1000; % sampling rate (Hz)
fP = logspace(log10(0.1),log10(10),21); % phase: frequencies (Hz)
fA = logspace(log10(10),log10(200),21); % amplitude: frequencies (Hz)

% load data
cd([ddir subj filesep tank]);
% dat = loadtxtfile(flnm)';
Nfl = length(dir('*_*.txt'));
dat = [];
for ii = 0:Nfl-1
    fnm = [tank(end-5:end-2) '_' num2str(ii) '.txt'];
    fid = fopen(fnm,'r');
    dat = [dat; fread(fid,inf,'uint16')];
    fclose(fid);
end
N = length(dat); % samples
t = (1:N)/dfs/3600; % hr
ind = find(t>=dhr(1) & t<=dhr(2));
dat = dat(ind); % limit to time window
    
% compute phase-amplitude modulation index (Tort et al 2010)
rP = [fP(1:end-1)' fP(2:end)'];
rA = [fA(1:end-1)' fA(2:end)'];
NP = size(rP,1);
NA = size(rA,1);
MI = zeros(NA,NP);
edges = -pi:delta:pi;
x = edges(1:end-1)+delta/2;
Nx = length(x); count = 0;
for iP = 1:NP
    [bP,aP] = butter(2,rP(iP,:)/(dfs/2));
    P = angle(hilbert(filtfilt(bP,aP,dat)));
    if cph
        ECDFt = sort(P);
        ECDFx = (1:length(ECDFt))/length(ECDFt);
        cx = interp1(ECDFt,ECDFx,P,'linear');
        P = 2*pi*cx-pi; % phase corrected for waveform asymmetry (Siapas et al 2005)
    end
    [~,bin] = histc(P,edges);
    for iA = 1:NA
        [bA,aA] = butter(2,rA(iA,:)/(dfs/2));
        A = abs(hilbert(filtfilt(bA,aA,dat)));
        Am = zeros(1,Nx);
        for jj = 1:Nx
            Am(jj) = mean(A(bin==jj));
        end
        Am = Am/sum(Am);
        Am(Am==0) = 1e-10; % to avoid -Inf for log(0)
        MI(iA,iP) = (log(Nx)+sum(Am.*log(Am)))/log(Nx);
        count = count + 1;
        disp([num2str(count) '/' num2str(NP*NA)]);
    end
end
    
% plot
figure('Name',tank,'NumberTitle','off','Units','normalized','Position',[1/4 1/4 1/2 1/2],'Color','w');
x = logspace(log10(fP(1)),log10(fP(end)),1000);
y = logspace(log10(fA(1)),log10(fA(end)),1000);
imagesc(1:length(x),1:length(y),smooth2(MI,2,2));
xtick = [0.1:0.1:1 2:10];
xticknm = {'0.1';'';'';'';'';'';'';'';'';'1';'';'';'';'';'';'';'';'';'10'};
xind = zeros(size(xtick));
for jj = 1:length(xtick)
    [~,xind(jj)] = min(abs(x-xtick(jj)));
end
ytick = [10:10:100 200];
yticknm = {'10';'';'';'';'';'';'';'';'';'100';''};
yind = zeros(size(ytick));
for jj = 1:length(ytick)
    [~,yind(jj)] = min(abs(y-ytick(jj)));
end
axis square; set(gca,'Box','off','YDir','normal','TickDir','out','XTick',xind,'XTickLabel',xticknm,'YTick',yind,'YTickLabel',yticknm);
xlabel('phase frequency (Hz)'); ylabel('amplitude frequency (Hz)');
title(tank,'Interpreter','none');
hc = colorbar; set(get(hc,'YLabel'),'String','modulation index');
