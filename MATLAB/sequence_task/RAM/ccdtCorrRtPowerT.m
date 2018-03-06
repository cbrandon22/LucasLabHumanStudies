function ccdtCorrRtPowerT
% function ccdtCorrRtPowerT
%   Neural correlates of CCDT: neural power in a chosen frequency band
%   across trials sorted by RT.
%
%   DR 03/2016

% parameters
ddir = 'F:\'; % directory
subj = 'fidel'; % subject
dday = [110415; 110615; 111115; 111315; 111615]; % day(s)
ch = 10; % channels
band = [5 10]; % frequency band (Hz)
sevt = 'cue'; % sort event: 'cue' or 'go'
cuewin = [-400 2500]; % cue-aligned window

% compile data across days
cd([ddir subj '\processed\ccdtCorr\']);
for ich = 1:length(ch)
    RTa = []; DTa = []; Pa = []; Pae = []; DTae = [];
    for iday = 1:length(dday)
        fls = dir([subj '_*_' num2str(dday(iday)) '_*_' num2str(ch(ich)) '.mat']);
        if length(fls)~=1, disp('file not found'); return; end
        load(fls(1).name);
        RTa = [RTa; RT]; DTa = [DTa; DT];
        Pa = [Pa; abs(Pcc)]; % amplitude
        DTae = [DTae; DTe];
        Pae = [Pae; abs(Pcce)]; % amplitude
    end
    
    % band
    Pa = squeeze(mean(Pa(:,freq>=band(1)&freq<=band(2),:),2)); % mean power in band
    Ptot = mean(Pa,2); iout = find(Ptot>prctile(Ptot,95)); % dismiss outlier trials (high power)
    Pa(iout,:) = []; RTa(iout) = []; DTa(iout) = [];
    Pae = squeeze(mean(Pae(:,freq>=band(1)&freq<=band(2),:),2));
    Ptot = mean(Pae,2); iout = find(Ptot>prctile(Ptot,95));
    Pae(iout,:) = []; DTae(iout) = [];
    
    % re-align
    if strcmp(sevt,'cue')
        if ich==1
            fs = 1000/mean(diff(t)); % samp/s
            swin = round(cuewin(1)/1000*fs):round(cuewin(2)/1000*fs);
        end
        Pc = zeros(size(Pa,1),length(swin));
        for ii = 1:size(Pa,1)
            [~,ind] = min(abs(t+DTa(ii)));
            Pc(ii,:) = Pa(ii,swin+ind);
        end
        Pa = Pc;
        t = swin/fs*1000;
    end
    
    % sort and smooth
    [RTa,isrt] = sort(RTa); DTa = DTa(isrt);
    Pa = Pa(isrt,:);
    Pa = smooth2(Pa,10,10);
    [DTae,isrte] = sort(DTae);
    Pae = Pae(isrte,:);
    Pae = smooth2(Pae,10,10);
    
    % plot
    figure('Name',num2str(ch(ich)),'NumberTitle','off','Units','normalized','Position',[1/8 1/6 3/4 2/3],'Color','w');
    subplot(1,2,1), hold on; colormap(parula(256));
    imagesc(t,1:length(RTa),Pa); hold on;
    axis_square(gca);
    set(gca,'Box','off','YDir','normal','TickDir','out','YLim',[1 length(RTa)],'XLim',[t(1) t(end)]);
    xlabel('time (s)'); ylabel('trials'); title({['channel ' num2str(ch(ich)) '  (' num2str(band(1)) '-' num2str(band(2)) 'Hz)'],[]});
    line([0 0],get(gca,'YLim'),'Color','w');
    text(0,max(get(gca,'YLim')),sevt,'HorizontalAlignment','center','VerticalAlignment','bottom');
    if strcmp(sevt,'go')
        plot(-DTa,1:length(RTa),'wo');
        text(-DTa(end),max(get(gca,'YLim')),'cue','HorizontalAlignment','center','VerticalAlignment','bottom');
        plot(RTa,1:length(RTa),'w+');
        text(RTa(end),max(get(gca,'YLim')),'move','HorizontalAlignment','center','VerticalAlignment','bottom');
    elseif strcmp(sevt,'cue')
        plot(DTa,1:length(RTa),'wo');
        text(DTa(end),max(get(gca,'YLim')),'go','HorizontalAlignment','center','VerticalAlignment','bottom');
        plot(DTa+RTa,1:length(RTa),'w+');
        text(DTa(end)+RTa(end),max(get(gca,'YLim')),'move','HorizontalAlignment','center','VerticalAlignment','bottom');
%         plot(RTa,1:length(RTa),'w--');
%         text(RTa(end),max(get(gca,'YLim')),'cue+RT','HorizontalAlignment','center','VerticalAlignment','bottom');
    end
    clim = get(gca,'CLim');
    pos = get(gca,'Position'); hcb = colorbar('East'); cbpos = get(hcb,'Position');
    set(hcb,'Position',[pos(1)+pos(3)+cbpos(3)/2, cbpos(2:4)],'YAxisLocation','right');
    ylabel(hcb,'amplitude');
    subplot(1,2,2), hold on; colormap(parula(256));
    imagesc(te,1:size(Pae,1),Pae); hold on;
    axis_square(gca);
    set(gca,'Box','off','YDir','normal','TickDir','out','YLim',[1 length(RTa)],'XLim',cuewin,'CLim',clim);
    xlabel('time (s)'); ylabel('trials'); title({'error trials',[]});
    line([0 0],get(gca,'YLim'),'Color','w');
    text(0,size(Pae,1),'cue','HorizontalAlignment','center','VerticalAlignment','bottom');
    plot(DTae,1:size(Pae,1),'w+');
    text(DTae(end),size(Pae,1),'move','HorizontalAlignment','center','VerticalAlignment','bottom');
    pos = get(gca,'Position'); hcb = colorbar('East'); cbpos = get(hcb,'Position');
    set(hcb,'Position',[pos(1)+pos(3)+cbpos(3)/2, cbpos(2:4)],'YAxisLocation','right');
    ylabel(hcb,'amplitude');
    if length(ch)>1
        orient landscape
        print('ccdtCorrRtPowerT.ps','-dpsc2','-append');
        close
    end
end
