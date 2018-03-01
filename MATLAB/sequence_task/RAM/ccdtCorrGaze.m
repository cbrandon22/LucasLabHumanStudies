function ccdtCorrGaze
% function ccdtCorrGaze
%   Neural correlates of CCDT: gaze distance from target.
%
%   DR 03/2016

% parameters
ddir = 'F:\'; % directory
subj = 'fidel'; % subject
dday = [110415; 110615; 111115; 111315; 111615]; % day(s)
acue = 1; % align to go (0) or cue (1)?
twin = [-300 1200]; % per-cue time window (ms)
pci = 1; % plot confidence intervals? (0 or 1)
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
GDa = []; RTa = []; DTa = [];
for iday = 1:length(dday)
    fls = dir([subj '_*_' num2str(dday(iday)) '_*_1.mat']); % all channels have gaze info
    if length(fls)~=1, disp('file not found'); return; end
    load(fls(1).name);
    RTa = [RTa; RT];
    DTa = [DTa; DT];
    
    % gaze distance
    GD = zeros(size(eyeXcc));
    for ii = 1:length(targ)
        itarg = targ(ii)-210;
        GD(ii,:) = sqrt((eyeXcc(ii,:)-Px(itarg)).^2+(eyeYcc(ii,:)-Py(itarg)).^2);
    end
    GDa = [GDa; GD];
end
GDa = (GDa-50)/ppi*2.54; % cm + hack

% re-align to cue
if acue
    fs = 1000/mean(diff(t)); % samp/s
    swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
    GDc = zeros(size(GDa,1),length(swin));
    for ii = 1:size(GDa,1)
        [~,ind] = min(abs(t+DTa(ii)));
        GDc(ii,:) = GDa(ii,swin+ind);
    end
    GDa = GDc;
    t = swin/fs*1000;
end

% correlation with RT
[C,Cp] = corr(GDa,RTa,'type','Spearman');

% bootstrap 95% CI on median distance
M = median(GDa);
if pci
    CI = zeros(2,length(M));
    for ii = 1:length(M)
        CI(:,ii) = bootci(1000,@(x)(median(x)),GDa(:,ii));
    end
end
    
% plot
figure('Name','ccdtCorrGaze','NumberTitle','off','Units','normalized','Position',[1/3 1/8 1/3 3/4],'Color','w');
subplot(2,1,1), hold on;
if pci, patch([t fliplr(t)],[CI(1,:) fliplr(CI(2,:))],'k','FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); end
plot(t,M,'k');
axis square; set(gca,'Box','off'); if acue, set(gca,'XLim',[t(1) t(end)]); xlabel('time from cue (ms)'); else set(gca,'XLim',[t(1) 0]); xlabel('time from go (ms)'); end
ylabel('gaze error (cm)'); title(['N = ' num2str(size(GDa,1))]);
if ~acue
    line(-1*[min(DT) min(DT)],get(gca,'YLim'),'Color','k','LineStyle',':');
    line(-1*[max(DT) max(DT)],get(gca,'YLim'),'Color','k','LineStyle',':');
    text(-1*mean(DT),max(get(gca,'YLim')),'cue (range)','HorizontalAlignment','center','VerticalAlignment','bottom');
else
    line([min(DT) min(DT)],get(gca,'YLim'),'Color','k','LineStyle',':');
    text(min(DT),max(get(gca,'YLim')),'go-->','HorizontalAlignment','left','VerticalAlignment','bottom');    
end
subplot(2,1,2), hold on;
plot(t,C,'k');
axis square; set(gca,'Box','off'); if acue, set(gca,'XLim',[t(1) t(end)]); xlabel('time from cue (ms)'); else set(gca,'XLim',[t(1) 0]); xlabel('time from go (ms)'); end
Cp(Cp>0.01) = NaN; Cp(~isnan(Cp)) = max(get(gca,'YLim'));
plot(t,Cp,'r','LineWidth',2);
ylabel('Spearman''s rank correlation coefficient'); title('correlation between gaze error and reaction time');
if ~acue
    line(-1*[min(DT) min(DT)],get(gca,'YLim'),'Color','k','LineStyle',':');
    line(-1*[max(DT) max(DT)],get(gca,'YLim'),'Color','k','LineStyle',':');
else
    line([min(DT) min(DT)],get(gca,'YLim'),'Color','k','LineStyle',':');
end
