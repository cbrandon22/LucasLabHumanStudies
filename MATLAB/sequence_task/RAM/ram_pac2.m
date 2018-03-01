function ram_pac2
% function ram_pac2
%   Compute phase-amplitude coupling from recording files (no stim):
%   comodulogram plot
% 
%   DR 03/2015

% parameters
ddir = 'G:\ditka\'; % directory
dtank = 'ditka_rec_DT1_032715'; % data tank
dblock = 1; % data block
ch = 3; % channel(s)
delta = pi/8; % phase bin size (rad)
cph = 1; % correct for waveform asymmetry? (0 or 1)

% constants
fs = 24414.06; % sampling rate (Hz)
fP = logspace(log10(0.1),log10(10),21); % phase: frequencies (Hz)
fA = logspace(log10(10),log10(300),21); % amplitude: frequencies (Hz)

% load data
cd([ddir dtank '\Block-' num2str(dblock)]);
clims = zeros(length(ch),2); hf = zeros(length(ch),1);
for ii = 1:length(ch)
    fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(ch(ii)) '.sev'],'r');
    fseek(fid,40,'bof');
    dat = fread(fid,inf,'single');
    fclose(fid);
    [b,a] = butter(2,500/(fs/2)); % antialiasing
    dat = filtfilt(b,a,dat);
    dat = downsample(dat',10)'; % downsample
    dfs = round(fs/10);
    
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
    hf(ii) = figure('Name',dtank,'NumberTitle','off','Units','normalized','Position',[1/4 1/4 1/2 1/2],'Color','w');
    x = logspace(log10(fP(1)),log10(fP(end)),1000);
    y = logspace(log10(fA(1)),log10(fA(end)),1000);
    imagesc(1:length(x),1:length(y),smooth2(MI,2,2)); clims(ii,:) = get(gca,'CLim');
    xtick = [0.1:0.1:1 2:10];
    xticknm = {'0.1';'';'';'';'';'';'';'';'';'1';'';'';'';'';'';'';'';'';'10'};
    xind = zeros(size(xtick));
    for jj = 1:length(xtick)
        [~,xind(jj)] = min(abs(x-xtick(jj)));
    end
    ytick = [10:10:100 200 300];
    yticknm = {'10';'';'';'';'';'';'';'';'';'100';'';''};
    yind = zeros(size(ytick));
    for jj = 1:length(ytick)
        [~,yind(jj)] = min(abs(y-ytick(jj)));
    end
    axis square; set(gca,'Box','off','YDir','normal','TickDir','out','XTick',xind,'XTickLabel',xticknm,'YTick',yind,'YTickLabel',yticknm);
    xlabel('phase frequency (Hz)'); ylabel('amplitude frequency (Hz)');
    title([dtank ', block ' num2str(dblock) ', ch ' num2str(ch(ii))],'Interpreter','none');
    hc = colorbar; set(get(hc,'YLabel'),'String','modulation index');
end
if length(ch)>1
    clmx = [min(clims(:,1)) max(clims(:,2))];
    for ii = 1:length(ch)
        figure(hf(ii));
        set(gca,'CLim',clmx);
        print([ddir dtank(end-5:end-2) '_pac2.ps'],'-dpsc2','-append');
    end
end