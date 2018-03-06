function ccdtCorrRtPower
% function ccdtCorrRtPower
%   Neural correlates of CCDT: neural power in a chosen frequency band and
%   time window as a function of RT.
%
%   DR 03/2016

% parameters
ddir = '/Users/tnl/Desktop/'; % directory
subj = 'fidel'; % subject
dday = [110415; 110615; 111115; 111315; 111615]; % day(s)
ch = 1:16; % channels
band = [5 12]; % frequency band (Hz)
twin = [-400 0]; % cue-aligned window (ms)

% compile data
cd([ddir subj '/processed/ccdtCorr/']);
CC = zeros(length(ch),2);
for ich = 1:length(ch)
    RTa = []; DTa = []; Pa = [];
    for iday = 1:length(dday)
        fls = dir([subj '_*_' num2str(dday(iday)) '_*_' num2str(ch(ich)) '.mat']);
        if length(fls)~=1, disp('file not found'); return; end
        load(fls(1).name);
        RTa = [RTa; RT]; DTa = [DTa; DT];
        Pa = [Pa; abs(Pcc)]; % amplitude
    end
    
    % window on cue event and compute average power
    if ich==1
        fs = 1000/mean(diff(t)); % samp/s
        swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
        ifreq = find(freq>=band(1)&freq<=band(2));
    end
    Pavg = zeros(size(Pa,1),1);
    for ii = 1:size(Pa,1)
        [~,ind] = min(abs(t+DTa(ii)));
        c = squeeze(Pa(ii,ifreq,swin+ind));
        Pavg(ii) = trimmean(c(:),5);
    end
    
    % correlation with RT
    [CC(ich,1),CC(ich,2)] = corr(Pavg,RTa,'type','Pearson');
end

% plot
figure('Name','ccdtCorrRTPower','NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
plot(ch,CC(:,1),'ko','MarkerFaceColor','k'); hold on;
for ich = 1:length(ch)
    if CC(ich,2)<0.001
        text(ch(ich),CC(ich,1),'*','HorizontalAlignment','center','VerticalAlignment','top','FontSize',14);
    end
end
axis square; set(gca,'Box','off','XTick',ch);
xlabel('channel'); ylabel('power-RT correlation'); title([num2str(band(1)) '-' num2str(band(2)) ' Hz, ' num2str(twin(1)) '-' num2str(twin(2)) ' ms from cue']);
