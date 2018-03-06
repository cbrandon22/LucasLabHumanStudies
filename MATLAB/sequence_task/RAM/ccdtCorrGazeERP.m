function ccdtCorrGazeERP
% function ccdtCorrGazeERP
%   Neural correlates of CCDT: event-related potentials at cue onset as a
%   function of gaze distance.
%
%   DR 03/2016

% parameters
ddir = 'F:\'; % directory
subj = 'fidel'; % subject
dday = [110415; 110615; 111115; 111315; 111615]; % day(s)
ch = 1:16; % channels
twin = [-200 500]; % time window relative to cue (ms)
lpf = 50; % low-pass filter cutoff (Hz, no filter if empty)
gp = [20 80]; % gaze distribution groupings (percentile)
ppi = 92; % pixels per inch

% constants (same as in VPLT.m)
im_pos = [727, 281, 535, 535]; % image position on monitor (pixels)
vt_size = [15 15]; % visual target size (pixels)
Px = [im_pos(1)-vt_size(1)/2; im_pos(1)+im_pos(3)+vt_size(1)/2; im_pos(1)-vt_size(1)/2;...
      im_pos(1)+im_pos(3)+vt_size(1)/2; im_pos(1)+im_pos(3)/2; im_pos(1)+im_pos(3)/4;...
      im_pos(1)+im_pos(3)/2; im_pos(1)+3*im_pos(3)/4; im_pos(1)+im_pos(3)/2]; % X pixel location of center of calibration targets
Py = [im_pos(2)-vt_size(2)/2; im_pos(2)-vt_size(2)/2; im_pos(2)+im_pos(4)+vt_size(1)/2;
      im_pos(2)+im_pos(4)+vt_size(1)/2; im_pos(2)+im_pos(4)/2; im_pos(2)+im_pos(4)/2;...
      im_pos(2)+im_pos(4)/4; im_pos(2)+im_pos(4)/2; im_pos(2)+3*im_pos(4)/4]; % Y pixel location of center of calibration targets

% compile data across days
cd([ddir subj '\processed\ccdtCorr\']);
for ich = 1:length(ch)
    RTa = []; DTa = []; Xa = []; GDa = [];
    for iday = 1:length(dday)
        fls = dir([subj '_*_' num2str(dday(iday)) '_*_' num2str(ch(ich)) '.mat']);
        if length(fls)~=1, disp('file not found'); return; end
        load(fls(1).name);
        GD = zeros(size(eyeXcc));
        for ii = 1:length(targ)
            itarg = targ(ii)-210;
            GD(ii,:) = sqrt((eyeXcc(ii,:)-Px(itarg)).^2+(eyeYcc(ii,:)-Py(itarg)).^2);
        end
        RTa = [RTa; RT]; DTa = [DTa; DT]; Xa = [Xa; Xcc]; GDa = [GDa; GD];
    end
    GDa = GDa/ppi*2.54; % cm
    if ich==1
        fs = 1000/mean(diff(t)); % samp/s
        swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
        if ~isempty(lpf), [b,a] = butter(2,lpf/(fs/2)); end
        GDac = zeros(size(GDa,1),1); % gaze distance from cue onset to cue+500ms
        for ii = 1:length(GDac)
            [~,ind] = min(abs(t+DTa(ii)));
            GDac(ii) = mean(GDa(ii,ind:ind+round(0.5*fs)));
        end
        [GDac,isrt] = sort(GDac);
    end
    if ~isempty(lpf), Xa = filtfilt(b,a,Xa')'; end
        
    % cue aligned, sorted by gaze distance
    Xc = zeros(size(Xa,1),length(swin));
    for ii = 1:size(Xa,1)
        [~,ind] = min(abs(t+DTa(ii)));
        Xc(ii,:) = Xa(ii,swin+ind);
    end
    Xc = Xc(isrt,:);
    ta = sum(abs(Xc),2); ibad = find(ta>prctile(ta,98)); Xc(ibad,:) = []; GDp = GDac; GDp(ibad) = []; RTa(ibad) = []; % remove outliers
    Nt = size(Xc,1);
        
    % plot
    t = swin/fs*1000;
    figure('Name','ccdtCorrGazeERP','NumberTitle','off','Units','normalized','Position',[1/6 1/6 2/3 2/3],'Color','w');
    subplot(3,3,1:2); hold on;
    i1 = find(GDp<prctile(GDp,gp(1)));
    plot(t,mean(Xc(i1,:)),'k','LineWidth',2);
    i2 = find(GDp>prctile(GDp,gp(2)));
    plot(t,mean(Xc(i2,:)),'Color',[0.8 0.8 0.8],'LineWidth',2);
    set(gca,'Box','off','XLim',[t(1) t(end)]);
    ylabel('\muV');
    subplot(3,3,3), hold on;
    plot(GDp,RTa,'k.'); axis tight;
    set(gca,'Box','off');
    cc = corr(GDp,RTa,'type','Spearman'); text(max(get(gca,'XLim')),max(get(gca,'YLim')),['cc = ' num2str(cc,'%4.2f')],'HorizontalAlignment','right','VerticalAlignment','top');
    ylabel('RT (ms)'); xlabel('gaze distance (cm)');
    subplot(3,3,[4 5 7 8]), hold on;
    imagesc(t,1:Nt,smooth2(Xc,10,10));
    set(gca,'Box','off','YDir','normal','TickDir','out','YLim',[1 Nt],'XLim',[t(1) t(end)]);
    line(t(end)*ones(1,2),[min(i1) max(i1)],'Color','k','LineWidth',4);
    line(t(end)*ones(1,2),[min(i2) max(i2)],'Color',[0.8 0.8 0.8],'LineWidth',4);
    xlabel('time from cue (ms)'); ylabel('trials'); title(num2str(ch(ich)));
    subplot(3,3,[6 9]), hold on;
    plot(GDp,1:Nt,'k','LineWidth',2); axis tight;
    set(gca,'Box','off','YLim',[1 Nt]);
    xlabel('gaze distance (cm)');
    if length(ch)>1
        orient landscape
        print('ccdtCorrGazeERP.ps','-dpsc2','-append');
        close
    end
end