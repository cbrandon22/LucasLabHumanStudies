function ccdtCorrRtERP
% function ccdtCorrRtERP
%   Neural correlates of CCDT: event-related potentials at cue onset as a
%   function of reaction time.
%
%   DR 03/2016

% parameters
ddir = 'F:\'; % directory
subj = 'fidel'; % subject
dday = [110415; 110615; 111115; 111315; 111615]; % day(s)
ch = 1:16; % channels
twin = [-200 500]; % time window relative to cue (ms)
lpf = 50; % low-pass filter cutoff (Hz, no filter if empty)
gp = [20 80]; % RT distribution groupings (percentile)

% compile data across days
cd([ddir subj '\processed\ccdtCorr\']);
for ich = 1:length(ch)
    RTa = []; DTa = []; Xa = [];
    for iday = 1:length(dday)
        fls = dir([subj '_*_' num2str(dday(iday)) '_*_' num2str(ch(ich)) '.mat']);
        if length(fls)~=1, disp('file not found'); return; end
        load(fls(1).name);
        RTa = [RTa; RT]; DTa = [DTa; DT]; Xa = [Xa; Xcc];
    end
    if ich==1
        fs = 1000/mean(diff(t)); % samp/s
        swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
        if ~isempty(lpf), [b,a] = butter(2,lpf/(fs/2)); end
    end
    if ~isempty(lpf), Xa = filtfilt(b,a,Xa')'; end
        
    % cue aligned, sorted by RT
    Xc = zeros(size(Xa,1),length(swin));
    for ii = 1:size(Xa,1)
        [~,ind] = min(abs(t+DTa(ii)));
        Xc(ii,:) = Xa(ii,swin+ind);
    end
    [RTa,isrt] = sort(RTa);
    Xc = Xc(isrt,:);
    ta = sum(abs(Xc),2); ibad = find(ta>prctile(ta,98)); Xc(ibad,:) = []; RTp = RTa; RTp(ibad) = []; % remove outliers
    Nt = size(Xc,1);
        
    % plot
    t = swin/fs*1000;
    figure('Name','ccdtCorrRtERP','NumberTitle','off','Units','normalized','Position',[1/6 1/6 2/3 2/3],'Color','w');
    subplot(3,3,1:2); hold on;
    i1 = find(RTp<prctile(RTp,gp(1)));
    plot(t,mean(Xc(i1,:)),'k','LineWidth',2);
    i2 = find(RTp>prctile(RTp,gp(2)));
    plot(t,mean(Xc(i2,:)),'Color',[0.8 0.8 0.8],'LineWidth',2);
    set(gca,'Box','off','XLim',[t(1) t(end)]);
    ylabel('\muV');
    subplot(3,3,[4 5 7 8]), hold on;
    imagesc(t,1:Nt,smooth2(Xc,10,10));
    set(gca,'Box','off','YDir','normal','TickDir','out','YLim',[1 Nt],'XLim',[t(1) t(end)]);
    line(t(end)*ones(1,2),[min(i1) max(i1)],'Color','k','LineWidth',4);
    line(t(end)*ones(1,2),[min(i2) max(i2)],'Color',[0.8 0.8 0.8],'LineWidth',4);
    xlabel('time from cue (ms)'); ylabel('trials'); title(num2str(ch(ich)));
    subplot(3,3,[6 9]), hold on;
    plot(RTp,1:Nt,'k','LineWidth',2); axis tight;
    set(gca,'Box','off','YLim',[1 Nt]);
    xlabel('RT (ms)');
    if length(ch)>1
        orient landscape
        print('ccdtCorrRtERP.ps','-dpsc2','-append');
        close
    end
end