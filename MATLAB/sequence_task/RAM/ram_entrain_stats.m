function ram_entrain_stats
% function ram_entrain_stats
%   Find amplitude and frequency of oscillatory potential evoked after each
%   stimulus train. Save to prepared data file.
%
%   DR 04/2015

% parameters
ddir = 'D:\'; % directory
subj = 'ditka'; % subject
tank = 'ditka_stim_DT1_033015'; % tank
block = 7; % block
ch = 12; % recording channel
twin = [-10 60]; % time window around train offset (ms)
lpf = 250; % lowpass filter cutoff to remove irrelevant local extrema
fullcyc = 1; % compute freq from full cycle (1) or half cycle (0) after reference point?
debug = 0; % debug plot? (0 or 1)

% window
try load([ddir subj '\processed\' tank '_' num2str(block) '.mat'],'-mat');
catch disp('run ''ram_stimtrig'' first'); return; end
if spuls>1, itrain = itrain + round((spuls-1)/sfreq*fs); end % last pulse of each train (sample)
swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
indwin = itrain*ones(1,length(swin)) + ones(length(itrain),1)*swin;
indwin(indwin(:,1)<1,:) = [];
N = size(indwin);
t = linspace(twin(1),twin(2),N(2));

% data
cd([ddir subj '\' tank '\Block-' num2str(block)]);
fid = fopen([tank '_Block-' num2str(block) '_xRAW_ch' num2str(ch) '.sev'],'r');
fseek(fid,40,'bof');
dat = fread(fid,inf,'single'); % uV
fclose(fid);
[b,a] = butter(2,lpf/(fs/2));
dat = filtfilt(b,a,dat); % lowpass filter
datwin = dat(indwin);
datwin = datwin - median(datwin,2)*ones(1,N(2)); % remove dc offset

% choose reference extremum from mean EP
figure
plot(mean(datwin),'k'); hold on;
axis tight; set(gca,'Box','off');
title('choose reference extremum');
[iref,~] = ginput(1);
close

% single trial EP stats
epstats = zeros(N(1),2); % [amplitude,frequency]
for jj = 1:N(1)
    cdat = datwin(jj,:);
    [~,imax,~,imin] = extrema(cdat);
    imax = sort(imax); imin = sort(imin);
    [smx,iimx] = min(abs(imax-iref)); % find extremum closest to reference extremum
    [smn,iimn] = min(abs(imin-iref));
    if smx < smn % reference extremum is a maximum
        ii1 = imax(iimx);
        iinext = find(imin > ii1);
        ii2 = imin(iinext(1)); % 1st minimum following the reference
        iinext = find(imax > ii1);
        ii3 = imax(iinext(1)); % 1st maximum following the reference 
    else % reference extremum is a minimum
        ii1 = imin(iimn);
        iinext = find(imax > ii1);
        ii2 = imax(iinext(1)); % 1st maximum following the reference
        iinext = find(imin > ii1);
        ii3 = imin(iinext(1)); % 1st minimum following the reference
    end
    if fullcyc
        epstats(jj,1) = abs(cdat(ii1)-cdat(ii2));
        epstats(jj,2) = fs/(ii3-ii1);
    else
        epstats(jj,1) = abs(cdat(ii1)-cdat(ii2));
        epstats(jj,2) = fs/(2*(ii2-ii1));
    end
    if debug
        figure
        plot(t,cdat,'k'); hold on;
        plot(t(ii1),cdat(ii1),'ro');
        plot(t(ii2),cdat(ii2),'bo');
        if fullcyc, plot(t(ii3),cdat(ii3),'go'); end
        set(gca,'Box','off','XLim',[t(1) t(end)]);
        title(['amp = ' num2str(epstats(jj,1),'%5.1f') ' \muV, freq = ' num2str(epstats(jj,2),'%4.1f') ' Hz']);
        pause; close;
    end
end

% save
entrain_stats{ch} = epstats;
save([ddir subj '\processed\' tank '_' num2str(block) '.mat'],'entrain_stats','-append','-mat');