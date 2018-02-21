function ccdtCorrRtPowerTF
% function ccdtCorrRtPowerTF
%   Neural correlates of CCDT: Comparison of neural spectrograms on short
%   and long RT trials.
%
%   DR 03/2016

% parameters
ddir = 'F:\'; % directory
subj = 'fidel'; % subject
dday = [110415; 110615; 111115; 111315; 111615]; % day(s)
ch = 4; % channels
rt1 = [0 250]; % RT group 1 (ms)
rt2 = [750 1000]; % RT group 2 (ms)
sevt = 'cue'; % sort event: 'cue' or 'go'
cuewin = [-400 2000]; % cue-aligned window
torp = 1; % plot tstat (0) or pvalue (1)?

% compile data across days
cd([ddir subj '\processed\ccdtCorr\']);
for ich = 1:length(ch)
    RTa = []; DTa = []; Pa = [];
    for iday = 1:length(dday)
        fls = dir([subj '_*_' num2str(dday(iday)) '_*_' num2str(ch(ich)) '.mat']);
        if length(fls)~=1, disp('file not found'); return; end
        load(fls(1).name);
        RTa = [RTa; RT]; DTa = [DTa; DT];
        Pa = [Pa; abs(Pcc)]; % LFP amplitude
    end
    
    % re-align
    if strcmp(sevt,'cue')
        if ich==1
            fs = 1000/mean(diff(t)); % samp/s
            swin = round(cuewin(1)/1000*fs):round(cuewin(2)/1000*fs);
        end
        Pc = zeros(size(Pa,1),size(Pa,2),length(swin));
        for ii = 1:size(Pa,1)
            [~,ind] = min(abs(t+DTa(ii)));
            Pc(ii,:,:) = Pa(ii,:,swin+ind);
        end
        Pa = Pc;
        t = swin/fs*1000;
    end
    
    % robust 2-sample t statistic
    i1 = find(RTa>=rt1(1) & RTa<=rt1(2));
    i2 = find(RTa>=rt2(1) & RTa<=rt2(2));
    Nf = length(freq); Nt = length(t);
    tstat = zeros(Nf,Nt);
    pstat = zeros(Nf,Nt);
    for ii = 1:Nf
        for jj = 1:Nt
            cP1 = squeeze(Pa(i1,ii,jj)); ibad1 = find(cP1>prctile(cP1,95)); cP1(ibad1) = []; N1 = length(cP1);
            cP2 = squeeze(Pa(i2,ii,jj)); ibad2 = find(cP2>prctile(cP2,95)); cP2(ibad2) = []; N2 = length(cP2);
            tstat(ii,jj) = (mean(cP1)-mean(cP2))/sqrt(var(cP1)/length(cP1)+var(cP2)/length(cP2));
            dof = (var(cP1)/N1+var(cP2)/N2)^2/((var(cP1)/N1)^2/(N1-1)+(var(cP2)/N2)^2/(N2-1));
            pstat(ii,jj) = 2*tcdf(-abs(tstat(ii,jj)),dof); % two-tailed p-value
        end
    end
    
    % plot
    if ich==1 % note: very specific to the 'freq' chosen in prep
        ytick = [1:10 20:10:100 200];
        yticknm = {'1';'';'';'';'';'';'';'';'';'10';'';'';'';'';'';'';'';'';'100';''};
        yind = zeros(size(ytick));
        for ii = 1:length(ytick)
            [~,yind(ii)] = min(abs(freq-ytick(ii)));
        end
    end
    figure('Name',num2str(ch(ich)),'NumberTitle','off','Units','normalized','Position',[1/4 1/6 1/2 2/3],'Color','w');
    cmap = colormap(parula(256)); colormap(flipud(cmap));
    tstat = smooth2(tstat,10,2); pstat = smooth2(pstat,10,2);
    if torp
        imagesc(t,1:length(freq),log10(pstat),[log10(0.01) 0]);
    else
        imagesc(t,1:length(freq),tstat,[-3 3]);
    end
    axis_square(gca);
    set(gca,'Box','off','YDir','normal','TickDir','out','YLim',[1 length(freq)],'YTick',yind,'YTickLabel',yticknm);
    xlabel('time (s)'); ylabel('frequency (Hz)'); title({['channel ' num2str(ch(ich))],[]});
    line([0 0],get(gca,'YLim'),'Color','w');
    if strcmp(sevt,'go')
        text(0,max(get(gca,'YLim')),'go','HorizontalAlignment','center','VerticalAlignment','bottom');
        line(-1*[min(DTa) min(DTa)],get(gca,'YLim'),'Color','w');
        line(-1*[max(DTa) max(DTa)],get(gca,'YLim'),'Color','w');
        text(-1*mean(DTa),max(get(gca,'YLim')),'cue (range)','HorizontalAlignment','center','VerticalAlignment','bottom');
    else
        text(0,max(get(gca,'YLim')),'cue','HorizontalAlignment','center','VerticalAlignment','bottom');
        line([min(DTa) min(DTa)],get(gca,'YLim'),'Color','w');
        line([max(DTa) max(DTa)],get(gca,'YLim'),'Color','w');
        text(mean(DTa),max(get(gca,'YLim')),'go (range)','HorizontalAlignment','center','VerticalAlignment','bottom');
    end
    pos = get(gca,'Position'); hcb = colorbar('East'); cbpos = get(hcb,'Position');
    set(hcb,'Position',[pos(1)+pos(3)+cbpos(3)/2, cbpos(2:4)],'YAxisLocation','right');
    if torp, ylabel(hcb,'log(p value)'); else ylabel(hcb,'t statistic'); end
    if length(ch)>1
        orient landscape
        print('ccdtCorrRtPowerTF.ps','-dpsc2','-append');
        close
    end
end
