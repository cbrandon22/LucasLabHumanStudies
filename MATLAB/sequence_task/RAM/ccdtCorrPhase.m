function ccdtCorrPhase
% function ccdtCorrPhase
%   Neural correlates of CCDT: event-related phase uniformity.
%
%   DR 03/2016

% parameters
ddir = 'F:\'; % directory
subj = 'fidel'; % subject
dday = [110615; 111115; 111315; 111615]; % day(s)
ch = 12; % channels
twin = [-200 400]; % time window relativel to even (ms)

% compile data
cd([ddir subj '\processed\ccdtCorr\']);
for ich = 1:length(ch)
    RTa = []; DTa = []; Pa = [];
    for iday = 1:length(dday)
        fls = dir([subj '_*_' num2str(dday(iday)) '_*_' num2str(ch(ich)) '.mat']);
        if length(fls)~=1, disp('file not found'); return; end
        load(fls(1).name);
        RTa = [RTa; RT]; DTa = [DTa; DT]; Pa = [Pa; angle(Pcc)]; % LFP phase
    end
    if ich==1
        fs = 1000/mean(diff(t)); % samp/s
        swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
    end
    Nt = length(RT); % trials
    Nf = length(freq); % frequencies
    Ns = length(swin); % samples
    
    % go-aligned phase
    [~,ind] = min(abs(t));
    cP = Pa(:,:,swin+ind);
    Pgo = zeros(Nf,Ns);
    for ii = 1:Nf
        for jj = 1:Ns
            Pgo(ii,jj) = rayleightest(squeeze(cP(:,ii,jj)));
        end
    end
    Pgo = log10(Pgo);
    Pgo = smooth2(Pgo,5,2);
    
    % cue-aligned phase
    cP = zeros(Nt,Nf,Ns);
    for ii = 1:Nt
        [~,ind] = min(abs(t+DTa(ii)));
        cP(ii,:,:) = Pa(ii,:,swin+ind);
    end
    Pcue = zeros(Nf,Ns);
    for ii = 1:Nf
        for jj = 1:Ns
            Pcue(ii,jj) = rayleightest(squeeze(cP(:,ii,jj)));
        end
    end
    Pcue = log10(Pcue);
    Pcue = smooth2(Pcue,5,2);
    
    % move-aligned phase
    cP = zeros(Nt,Nf,Ns);
    for ii = 1:Nt
        [~,ind] = min(abs(t+RTa(ii)));
        cP(ii,:,:) = Pa(ii,:,swin+ind);
    end
    Pmove = zeros(Nf,Ns);
    for ii = 1:Nf
        for jj = 1:Ns
            Pmove(ii,jj) = rayleightest(squeeze(cP(:,ii,jj)));
        end
    end
    Pmove = log10(Pmove);
    Pmove = smooth2(Pmove,5,2);
    
    % plot
    if ich==1 % note: very specific to the 'freq' chosen in prep
        ytick = [1:10 20:10:100 200];
        yticknm = {'1';'';'';'';'';'';'';'';'';'10';'';'';'';'';'';'';'';'';'100';''};
        yind = zeros(size(ytick));
        for ii = 1:length(ytick)
            [~,yind(ii)] = min(abs(freq-ytick(ii)));
        end
    end
    t = swin/fs*1000;
    figure('Name',num2str(ch(ich)),'NumberTitle','off','Units','normalized','Position',[1/8 1/4 3/4 1/2],'Color','w');
    cmap = colormap(parula(256)); colormap(flipud(cmap));
    subplot(1,3,1), hold on;
    imagesc(t,1:length(freq),Pcue,[-2 0])
    axis square; set(gca,'Box','off','XLim',[t(1) t(end)],'YDir','normal','TickDir','out','YLim',[1 length(freq)],'YTick',yind,'YTickLabel',yticknm);
    xlabel('time from cue (ms)'); ylabel('frequency (Hz)');
    subplot(1,3,2), hold on;
    imagesc(t,1:length(freq),Pgo,[-2 0])
    axis square; set(gca,'Box','off','XLim',[t(1) t(end)],'YDir','normal','TickDir','out','YLim',[1 length(freq)],'YTick',yind,'YTickLabel',yticknm);
    xlabel('time from go (ms)'); title(num2str(ch(ich)));
    subplot(1,3,3), hold on;
    imagesc(t,1:length(freq),Pmove,[-2 0])
    axis square; set(gca,'Box','off','XLim',[t(1) t(end)],'YDir','normal','TickDir','out','YLim',[1 length(freq)],'YTick',yind,'YTickLabel',yticknm);
    xlabel('time from move (ms)');
    pos = get(gca,'Position'); hcb = colorbar('East'); cbpos = get(hcb,'Position');
    set(hcb,'Position',[pos(1)+pos(3)+cbpos(3)/2, cbpos(2:4)],'YAxisLocation','right');
    ylabel(hcb,'log(p value)');
    if length(ch)>1
        orient landscape
        print('ccdtCorrPower.ps','-dpsc2','-append');
        close
    end
end


