function ccdtCorrERP
% function ccdtCorrERP
%   Neural correlates of CCDT: event-related potentials.
%
%   DR 03/2016

% parameters
ddir = 'F:\'; % directory
subj = 'fidel'; % subject
dday = [110415; 110615; 111115; 111315; 111615]; % day(s)
ch = 1:16; % channels
twin = [-200 400]; % time window relativel to even (ms)
lpf = 50; % low-pass filter cutoff (Hz, no filter if empty)

% compile data
cd([ddir subj '\processed\ccdtCorr\']);
erpGO = []; erpCUE = []; erpMOV = [];
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
    
    % go event
    [~,ind] = min(abs(t));
    erpGO = [erpGO; trimmean(Xa(:,swin+ind),5)];
    
    % cue event
    Xc = zeros(size(Xa,1),length(swin));
    for ii = 1:size(Xa,1)
        [~,ind] = min(abs(t+DTa(ii)));
        Xc(ii,:) = Xa(ii,swin+ind);
    end
    erpCUE = [erpCUE; trimmean(Xc,5)];
    
    % move event
    Xm = zeros(size(Xa,1),length(swin));
    for ii = 1:size(Xa,1)
        [~,ind] = min(abs(t-RTa(ii)));
        Xm(ii,:) = Xa(ii,swin+ind);
    end
    erpMOV = [erpMOV; trimmean(Xm,5)];
end

% plot
ylim = [min([erpGO(:);erpCUE(:);erpMOV(:)]) max([erpGO(:);erpCUE(:);erpMOV(:)])];
t = swin/fs*1000;
figure('Name','ccdtCorrERP','NumberTitle','off','Units','normalized','Position',[1/8 1/4 3/4 1/2],'Color','w');
cmap = colormap(jet(length(ch)));
subplot(1,3,1), hold on;
for ich = 1:length(ch)
    plot(t,erpCUE(ich,:),'Color',cmap(ich,:));
end
axis square; set(gca,'Box','off','XLim',[t(1) t(end)],'YLim',ylim);
xlabel('time from cue (ms)'); ylabel('\muV'); title(['N = ' num2str(size(Xc,1))]);
subplot(1,3,2), hold on;
for ich = 1:length(ch)
    plot(t,erpGO(ich,:),'Color',cmap(ich,:));
end
axis square; set(gca,'Box','off','XLim',[t(1) t(end)],'YLim',ylim);
xlabel('time from go (ms)');
subplot(1,3,3), hold on;
for ich = 1:length(ch)
    plot(t,erpMOV(ich,:),'Color',cmap(ich,:));
end
axis square; set(gca,'Box','off','XLim',[t(1) t(end)],'YLim',ylim);
xlabel('time from move (ms)');
legend(num2str(ch'),'Location','east');
