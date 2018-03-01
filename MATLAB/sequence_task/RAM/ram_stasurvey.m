function ram_stasurvey
% function ram_stasurvey
%   Stimulus-triggered average across trials: survey of evoked potentials
%   on each of the 8 electrodes on one array due to stimulation between
%   each of the 7 adjacent pairs on another array.
%
%   DR 03/2015

% parameters
ddir = 'D:\ditka\'; % data directory
dtkbk = {'ditka_stim_DT1_032315',10;...  % stim 1-2 (data tank, block)
         'ditka_stim_DT1_032315',9;...  % stim 2-3
         'ditka_stim_DT1_032315',7;...  % stim 3-4
         'ditka_stim_DT1_032315',6;...  % stim 4-5
         'ditka_stim_DT1_032315',5;...  % stim 5-6
         'ditka_stim_DT1_032315',3;...  % stim 6-7
         'ditka_stim_DT1_032315',4};    % stim 7-8
chrec = 1:8; % rec channels (should be 8 channels)
twin = [-100 300]; % sta time window (ms)
scale = 200; % uV
pci = 0; % plot 95% confidence intervals? (0 or 1)

% main
fs = 24414.06; % sampling rate (Hz)
swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
M = zeros(7,8,length(swin)); V = zeros(7,8,length(swin));
if length(chrec)~=8, error('must specify exactly 8 recording channels'); end
if size(dtkbk,1)~=7, error('must specify exactly 7 stimulating blocks'); end
for jj = 1:7
    disp([num2str(jj) '/' num2str(size(dtkbk,1))]);
    
    % stimulus times
    cd([ddir dtkbk{jj,1} '\Block-' num2str(dtkbk{jj,2})]);
    fid = fopen([dtkbk{jj,1} '_Block-' num2str(dtkbk{jj,2}) '_xRAW_ch25.sev'],'r');
    fseek(fid,40,'bof');
    dat = fread(fid,inf,'single');
    fclose(fid);
    itrain = find(dat(1:end-1)<0.5 & dat(2:end)>0.5); % trigger for each pulse/burst
    itrain = itrain + round(0.3/1000*fs); % adjust trigger time for ~0.3-ms delay until stim voltage starts
    indwin = itrain*ones(1,length(swin)) + ones(length(itrain),1)*swin;
    indwin(indwin(:,1)<1,:) = [];
    indwin(indwin(:,end)>length(dat),:) = [];
    
    % rec channels
    for ii = 1:8
        fid = fopen([dtkbk{jj,1} '_Block-' num2str(dtkbk{jj,2}) '_xRAW_ch' num2str(chrec(ii)) '.sev'],'r');
        fseek(fid,40,'bof');
        dat = fread(fid,inf,'single'); % uV
        fclose(fid);
        datwin = dat(indwin);
        datwin = datwin - median(datwin,2)*ones(1,length(swin)); % remove dc offset
        M(jj,ii,:) = mean(datwin);
        V(jj,ii,:) = 1.96*std(datwin)/sqrt(size(datwin,1));
    end
end

% plot
figure('Name','sta survey','NumberTitle','off','Units','normalized','Position',[1/8 1/8 3/4 3/4],'Color','w');
t = linspace(-0.4,0.4,length(swin));
for jj = 1:7
    for ii = 1:8
        m = squeeze(M(jj,ii,:))'/scale+jj;
        v = squeeze(V(jj,ii,:))'/scale;
        if pci, patch([t+ii,fliplr(t+ii)],[m+v,fliplr(m-v)],'k','FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); hold on; end
        plot(t+ii,m,'k','LineWidth',2); hold on;
    end
end
axis tight; set(gca,'Box','off','XLim',[0.5 8.5],'XTick',1:8,'YLim',[0 8],'YTick',1:7,'YTickLabel',{'1-2','2-3','3-4','4-5','5-6','6-7','7-8'});
line(max(get(gca,'XLim'))*ones(1,2),[7 8],'Color','k');
text(max(get(gca,'XLim')),7.5,[num2str(scale) ' \muV'],'HorizontalAlignment','left','VerticalAlignment','middle');
line([max(get(gca,'XLim'))-(t(end)-t(1)) max(get(gca,'XLim'))],[8 8],'Color','k');
text(max(get(gca,'XLim'))-(t(end)-t(1))/2,8,[num2str(diff(twin)) ' ms'],'HorizontalAlignment','center','VerticalAlignment','bottom');
xlabel('recording channels'); ylabel('stimulating channels');