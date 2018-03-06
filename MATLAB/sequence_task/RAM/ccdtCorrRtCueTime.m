function ccdtCorrRTCueTime
% function ccdtCorrRTCueTime
%   Neural correlates of CCDT: time series centered on cue onset for slow
%   and fast RT trials.
%
%   DR 03/2016

% parameters
ddir = '/Users/tnl/Desktop/'; % directory
subj = 'fidel'; % subject
dday = [110415; 110615; 111115; 111315; 111615]; % day(s)
ch = 10; % channels
band = [1 25]; % frequency band (Hz; don't filter if empty)
Nt = 20; % number of lowest and highest RT trials to plot
twin = [-400 500]; % time window relativel to cue onset (ms)
scale = 250; % amplitude scale factor

% compile data across days
cd([ddir subj '/processed/ccdtCorr/']);
for ich = 1:length(ch)
    RTa = []; DTa = []; Xa = [];
    for iday = 1:length(dday)
        fls = dir([subj '_*_' num2str(dday(iday)) '_*_' num2str(ch(ich)) '.mat']);
        if length(fls)~=1, disp('file not found'); return; end
        load(fls(1).name);
        RTa = [RTa; RT]; DTa = [DTa; DT]; Xa = [Xa; Xcc];
    end
    
    % sort by RT, filter, re-align on cue, divide into low/high RT groups
    [RTa,isrt] = sort(RTa);
    Xa = Xa(isrt,:); DTa = DTa(isrt);
    if ich==1
        fs = 1000/mean(diff(t)); % samp/s
        swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
        if ~isempty(band)
            [b,a] = butter(2,band/(fs/2));
        end
    end
    if ~isempty(band)
        Xa = filtfilt(b,a,Xa')';
    end
    Xc = zeros(size(Xa,1),length(swin));
    for ii = 1:size(Xa,1)
        [~,ind] = min(abs(t+DTa(ii)));
        Xc(ii,:) = Xa(ii,swin+ind);
    end
    X1 = Xc(1:Nt,:); RT1 = mean(RTa(1:Nt));
    X2 = Xc(end-Nt+1:end,:); RT2 = mean(RTa(end-Nt+1:end));
    t = swin/fs*1000;
    
    % plot
    figure('Name',num2str(ch(ich)),'NumberTitle','off','Units','normalized','Position',[1/4 1/6 1/2 2/3],'Color','w');
    subplot(1,2,1), hold on;
    imagesc(t,1:Nt,smooth2(abs(X1),40,3)); clim1 = get(gca,'CLim');
    for ii = 1:Nt
        plot(t,X1(ii,:)/scale+ii,'k');
    end
    set(gca,'Box','off','XLim',[t(1) t(end)],'YLim',[0.5 Nt+0.5]);
    xlabel('time from cue (ms)'); ylabel('trials'); title(['mean RT = ' num2str(round(RT1))]);
    subplot(1,2,2), hold on;
    imagesc(t,1:Nt,smooth2(abs(X2),20,3)); clim2 = get(gca,'CLim');
    for ii = 1:Nt
        plot(t,X2(ii,:)/scale+ii,'k');
    end
    set(gca,'Box','off','XLim',[t(1) t(end)],'YLim',[0.5 Nt+0.5]);
    xlabel('time from cue (ms)'); title(['mean RT = ' num2str(round(RT2))]);
    clim = [min(clim1(1),clim2(1)) max(clim1(2),clim2(2))];
    subplot(1,2,1), set(gca,'CLim',clim);
    subplot(1,2,2), set(gca,'CLim',clim);
    if length(ch)>1
        orient landscape
        print('ccdtCorrSpectrogram.ps','-dpsc2','-append');
        close
    end
end
