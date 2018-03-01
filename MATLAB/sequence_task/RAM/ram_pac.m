function [MI,PV,PP] = ram_pac(varargin)
% function ram_pac
%   Compute phase-amplitude coupling from recording files (no stim).
% 
%   DR 03/2015

% parameters
ddir = 'G:\'; % directory
subj = 'ditka'; % subject
tank = 'ditka_rec_DT1_032715'; % data tank
block = 1; % data block
ch = 1:24; % channel(s)
fP = [0.5 2]; % bandpass cutoffs - phase range
fA = [35 45]; % bandpass cutoffs - amplitude range
delta = pi/8; % phase bin size (rad)
cph = 1; % correct for waveform asymmetry? (0 or 1)
if nargin
    v2struct(varargin{1});
end

% data
MI = zeros(1,length(ch));
PV = zeros(1,length(ch));
PP = zeros(1,length(ch));
for jj = 1:length(ch)
    fs = 24414.06; % sampling rate (Hz)
    cd([ddir subj '\' tank '\Block-' num2str(block)]);
    fid = fopen([tank '_Block-' num2str(block) '_xRAW_ch' num2str(ch(jj)) '.sev'],'r');
    fseek(fid,40,'bof');
    dat = fread(fid,inf,'single');
    fclose(fid);
    [b,a] = butter(2,500/(fs/2)); % antialiasing
    dat = filtfilt(b,a,dat);
    dat = downsample(dat',10)'; % downsample
    fs = round(fs/10);

    % phases and amplitudes
    [bP,aP] = butter(2,fP/(fs/2));
    an = hilbert(filtfilt(bP,aP,dat));
    P = angle(an);
    if cph
        ECDFt = sort(P);
        ECDFx = (1:length(ECDFt))/length(ECDFt);
        cx = interp1(ECDFt,ECDFx,P,'linear');
        P = 2*pi*cx-pi; % phase corrected for waveform asymmetry (Siapas et al 2005)
    end
    [bA,aA] = butter(2,fA/(fs/2));
    A = abs(hilbert(filtfilt(bA,aA,dat)));
%     A = A(randi(length(A),size(A))); % shuffle to evaluate significance

    % phase-binned amplitude
    edges = -pi:delta:pi;
    [~,bin] = histc(P,edges);
    x = edges(1:end-1)+delta/2;
    N = length(x);
    Am = zeros(1,N);
    for ii = 1:N
        Am(ii) = mean(A(bin==ii));
    end
        
    % modulation index (Tort et al 2010)
    Am = Am/sum(Am);
    MI(jj) = (log(N)+sum(Am.*log(Am)))/log(N);
    
    % significance of unimodal tuning (on mean curve??)
    [~,PV(jj)] = circsig([Am' x'],0.05,0,1000);
    
    % preferred phase
    c = nansum(Am.*cos(x));
    s = nansum(Am.*sin(x));
    PP(jj) = atan2(s,c);

    % plot
    figure('Name',tank,'NumberTitle','off','Units','normalized','Position',[1/4 1/4 1/2 1/2],'Color','w');
    h = bar(x*180/pi,Am,'FaceColor',[0 0 0],'EdgeColor','k','BarWidth',1.0);
    axis tight; axis square; set(gca,'Box','off'); ylabel([num2str(fA(1)) '-' num2str(fA(2)) ' Hz amplitude (normalized)']); xlabel([num2str(fP(1)) '-' num2str(fP(2)) ' Hz phase (deg)']);
    tit = ['MI = ' num2str(MI(jj),'%1.1e') ', p = ' num2str(PV(jj),'%1.1e')]; if PV(jj)<0.05, tit = [tit ', ' num2str(round(PP(jj)*180/pi)) ' deg']; end
    title({[tank ', block ' num2str(block) ', ch ' num2str(ch(jj))],tit},'VerticalAlignment','bottom','Interpreter','none');
    
    % print
    if length(ch)>1
        print([ddir tank(end-5:end-2) '_pac.ps'],'-dpsc2','-append');
        close
    end
end
if length(ch)>1
    figure('Name',tank,'NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
    plot(ch,MI,'ko','LineStyle','-','MarkerFaceColor','k');
    axis square; set(gca,'Box','off'); xlabel('channel'); ylabel('modulation index');
    print([ddir tank(end-5:end-2) '_pac.ps'],'-dpsc2','-append');
    figure('Name',tank,'NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
    theta = 0:pi/64:2*pi;
    plot(cos(theta),sin(theta),'k');
    for jj = 1:length(ch)
        cpp = PP(jj);
        if MI(jj)>5e-4
            line([0 cos(cpp)],[0 sin(cpp)],'Color','k','LineWidth',2); hold on;
            text(cos(cpp),sin(cpp),num2str(ch(jj)),'HorizontalAlignment','center','VerticalAlignment','bottom');
        end
    end
    axis equal; axis square; set(gca,'Visible','off');
    print([ddir tank(end-5:end-2) '_pac.ps'],'-dpsc2','-append');
end