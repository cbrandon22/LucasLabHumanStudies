function ram_freq_anesthesia2L
% 
%   Compute power spectrum from data.
% 
%  EH 10/2016

for el = [2]
% parameters
% [2:13 15:32 34:36 38:40 42:46 48:54 56:64]; % Grid
ddir = 'L:\C\data'; % directory 'E:\data\C\data\eeg\' 'L:\C\data\'
subj = 'HUP119_i'; % subject 'HUP117_i'
dn = 0; % denoise the data (value is SD to use as cutoff), leave blank to not denoise data
elec = el; % Electrode number
dobehave = 1; % process behavior if = 1
gb = 1; % plot behavior (graph behavior)
snipit = 0; % time, in minutes to clip at the beginning of the spectrogram

tank = num2str(elec); % Electrode number 
% fs = 1000; % sampling rate (Hz)
nch = 1; % number of channels

% etime = 10; % length of time in minutes
% cor = [3.92 15.233]; % correct range
% LWM = 10.519; % LWM (loss of memory) in minutes from events file (lfpoffset / samplerate / 60) samplerate from params.txt
% LOC = 15.756; % LOC in minutes from events file
cd([ddir filesep subj filesep 'lfp.noreref']);
fname = ls('*.001');
[~,subjF,~] = fileparts(fname); % automatically obtain the correct filename for electrodes


sr = eegparams('samplerate',[ddir  filesep subj filesep 'lfp.noreref' filesep 'params.txt']); %obtain sample rate from params text file
timefix = sr*60; % convert sample rate to samples/min
% fs = sr;

events = [];
load([ddir filesep subj filesep 'behavioral' filesep 'session_0' filesep 'events.mat']);
rootF = [ddir filesep subj filesep 'lfp.noreref' filesep];
% etime = events(length(events)-1).lfpoffset/timefix; %obtain length of behavioral event trials (last trial is usually incomplete, hence length(events)-1
% etime = (events(length(events)-1).lfpoffset-events(1).lfpoffset)/timefix;
% etime = 11;
% subj = 'HUP119_i_24Jun16_0959'

for n=1:length(events)
    events(n).lfpfile = [rootF subjF];
end

if dobehave == 1
for n=1:length(events)
    events(n).lfpoffset2 = events(n).lfpoffset/timefix;
    events(n).lfpoffsetS = events(n).lfpoffset/sr;
end

    
% [correctR incorrectR] = eventStrategize(events); % find correct, incorrect, and pass (LOC) recall trials
% cor = [correctR(1).lfpoffset/timefix correctR(size(correctR,2)).lfpoffset/timefix];
% LWM = incorrectR(1).lfpoffset/timefix;
% LOC = correctR(size(correctR,2)).lfpoffset/timefix; % estmated based on final correct trial

ci = 1;
ii = 1;
rti = 1;
rti1 = 1;
resorte = 0; %if you're taking raw data, use 0 for the start of the ieeg data, otherwise use events(1).lfpoffset2 to resort if extracting only event-related data
yRT1 = []; xpc = []; ypc = []; xpcR = []; ypcR = []; yRT = []; yRT2 = []; xRT1 = []; yRT3 = [];

for i = 1:length(events)
    if i+2<length(events) && strcmp(events(i).type,'RECALL') && strcmp(events(i+2).type,'RECALL');
        events(i).correct = strcmp(events(i).response,events(i).target);
        events(i+1).correct = strcmp(events(i+1).response,events(i+1).target);
        events(i+2).correct = strcmp(events(i+2).response,events(i+2).target);
        events(i).percentcorrect = (events(i).correct + events(i+2).correct + events(i+1).correct)/3*100;
        events(i+1).percentcorrect = events(i).percentcorrect;
        events(i+2).percentcorrect = events(i).percentcorrect;
        xpc(ci) = events(i).lfpoffset2-resorte;
        xpc(ci+1) = events(i+1).lfpoffset2-resorte;
        xpc(ci+2) = events(i+2).lfpoffset2-resorte;
        ypc(ci) = events(i).percentcorrect;
        ypc(ci+1) = events(i+1).percentcorrect;
        ypc(ci+2) = events(i+2).percentcorrect;
        events(i).RT1 = events(i+1).lfpoffset-events(i).lfpoffset;
        events(i).RT2 = events(i+2).lfpoffset-events(i+1).lfpoffset;
        xRT1(rti1) = events(i).lfpoffset2-resorte;
        yRT1(rti1) = mean([events(i).RT1 events(i).RT2])/sr;
        rti1 = rti1+1;
        ci = ci + 3;   
    end
    if i+4<length(events) && strcmp(events(i).type,'RESP') && strcmp(events(i+2).type,'RESP') && strcmp(events(i+4).type,'RESP');
        events(i).correctRESP = strcmp(events(i).response,events(i).target);
        events(i+2).correctRESP = strcmp(events(i+2).response,events(i+2).target);
        events(i+4).correctRESP = strcmp(events(i+4).response,events(i+4).target);
        events(i).percentcorrectRESP = (events(i).correctRESP + events(i+2).correctRESP + events(i+4).correctRESP)/3*100;
        events(i+2).percentcorrectRESP = events(i).percentcorrectRESP;
        events(i+4).percentcorrectRESP = events(i).percentcorrectRESP;
        xpcR(ii) = events(i).lfpoffset2-resorte;
        xpcR(ii+1) = events(i+2).lfpoffset2-resorte;
        xpcR(ii+2) = events(i+4).lfpoffset2-resorte;
        ypcR(ii) = events(i).percentcorrectRESP;
        ypcR(ii+1) = events(i+2).percentcorrectRESP;
        ypcR(ii+2) = events(i+4).percentcorrectRESP;
        ii = ii + 3;   
    end
        if i+1<length(events) && strcmp(events(i).type,'SOUND') && strcmp(events(i+1).type,'RESP')
            events(i).RT = events(i+1).lfpoffset-events(i).lfpoffset;
            xRT(rti) = events(i).lfpoffset2-resorte;
            yRT(rti) = events(i).RT/sr;
            rti = rti + 1;
        
        end
        
end

    yRT2 = yRT/max(yRT)*100; % percent maximal RT
    if yRT1
    yRT3 = yRT1/max(yRT1)*100; % percent maximal RT
    end
end
% events = events(1); % take only the first event to look at entire time sequence.
% events.lfpoffset = 0;
% events.lfpfile = [ddir subj filesep 'lfp.noreref' filesep subjF];


% data
freq = logspace(log10(0.1),log10(200),100)'; % frequencies (Hz)
winsz = 5; % window size for spectral estimates (s) - determines frequency resolution
% try load([ddir subj filesep 'processed' filesep subj '_' tank '.mat'],'-mat');
% catch
%     DAT = [];   etime*60*1000
    
    if size(elec,2) == 1 % not bipolar, proceed normally
        channelSuffix = sprintf('.%03i',elec);
        [DAT,fs] = loadeData([ddir filesep subj filesep 'lfp.noreref' filesep subjF channelSuffix]);
    % remove noise (> or < # zscore SD from mean data)
    if dn > 0
        zDat = zscore(DAT);
        zR = find(zDat>dn | zDat < -dn);
        DAT(zR) = NaN;
    end
    elseif size(elec,2) == 2 % then it is bipolar
        channelSuffix1 = sprintf('.%03i',elec(1));
        channelSuffix2 = sprintf('.%03i',elec(2));
        [DAT1,fs1] = loadeData([ddir filesep subj filesep 'lfp.noreref' filesep subjF channelSuffix1]);
        [DAT2,fs2] = loadeData([ddir filesep subj filesep 'lfp.noreref' filesep subjF channelSuffix2]);
        DAT = DAT1 - DAT2;
        fs = fs1;
        % remove noise (> or < # zscore SD from mean data)
        if dn > 0
            zDat = zscore(DAT);
            zR = find(zDat>dn | zDat < -dn);
            DAT(zR) = NaN;
        end
    end
    
%     [DAT,rDAT] = an_getlfp_ms_wrapper(elec,events,etime*60*1000,0,0,[0 200],'bandpass',1); % resampled to 1000Hz
%     DAT = gete(elec,events,0);
%     DAT = geteeg_ms(elec,events,etime*60*1000,0,0);
%     DAT = an_getlfp_ms_wrapper(elec,events,etime*60*1000,0,0);

    
    disp('loading data...');
 
    
%     lastfullcycle = floor(length(DAT)/nch)*nch;
%     DAT(lastfullcycle+1:end) = [];
    DAT = reshape(DAT,nch,[])';
    N = size(DAT,1); % samples
    S = fix(N/fs)-(winsz-1); % seconds
    
    % spectra
    PS = zeros(nch,S,length(freq));
    fprintf('computing spectra...');
    for ii = 1:S % compute in 'winsz'-sec windows in 1-sec steps
        fprintf('%02d%%',fix(ii/S*100));
        ind = round((ii-1)*fs+1:(ii+winsz-1)*fs);
        for jj = 1:nch
            PS(jj,ii,:) = periodogram(DAT(ind,jj),hamming(length(ind)),freq,fs); % TODO: remove mean before spectral analysis?
        end
        if ii~=S, fprintf('\b\b\b'); end
    end
    fprintf('\n');
    
    % smoothing
    M = fix(S/60);
    PMr = zeros(nch,size(PS,3),M);
    PM = zeros(nch,size(PS,3),M);
    disp('smoothing spectra...');
    for ii = 1:M
        ind = (ii-1)*60+1:ii*60;
        for jj = 1:nch
            PMr(jj,:,ii) = trimmean(squeeze(PS(jj,ind,:)),10); % robust mean for every 60-sec (removes upper and lower 'winsz' seconds of each minute of data)
        end
    end
    x = 1:M;


    

    
    % end

% pre-plot
ind = find(x>snipit);
x = x(ind);
PM = PMr(:,:,ind);
for ich = 1:size(PM,1)
    PM(ich,:,:) = smooth2(squeeze(PM(ich,:,:)),2,1);
    PM(ich,:,:) = zscore(squeeze(PM(ich,:,:))')'; % zscore
end

    % save
    cd_mkdir([ddir filesep subj filesep 'processed']);
    save([ddir filesep subj filesep 'processed' filesep 'LGrid' filesep subj '_' tank],'PM','PMr','x','freq','DAT','PS','sr');
    if dobehave == 1
    save([ddir filesep subj filesep 'processed' filesep subj '_' tank 'behave'],'xpc','ypc','xpcR','ypcR','xRT','yRT','yRT2','xRT1','yRT1','yRT3');
    end

% plot
for ich = 1:size(PM,1)
    figure('Name',[subj ' Electrode: ' tank],'NumberTitle','off','Units','normalized','Position',[1/8 1/4 3/4 1/2],'Color','w');
    imagesc(x,1:length(freq),squeeze(PM(ich,:,:)));
    ytick = [1:10 20:10:100 200];
    yticknm = {'1';'';'';'';'';'';'';'';'';'10';'';'';'';'';'';'';'';'';'100';''};
    yind = zeros(size(ytick));
    for ii = 1:length(ytick)
        [~,yind(ii)] = min(abs(freq-ytick(ii)));
    end
    set(gca,'Box','off','YDir','normal','TickDir','out','YLim',[1 length(freq)],'XLim',[x(1) x(end)],'YTick',yind,'YTickLabel',yticknm);
    xlabel('time (min)'); ylabel('frequency (Hz)'); title([subj ' Electrode: ' tank],'Interpreter','none');
    set(gca,'CLim',[-3 3]);
    hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
%     line([LWM LWM],[0 200],'LineWidth',4,'Color','k')
%     line([LOC LOC],[0 200],'LineWidth',4,'Color','r')
%     line(cor,[1 1],'LineWidth',2,'Color','r')
    
if dobehave == 1 && gb == 1
hold on 
plot(xpc,ypc,'--r','LineWidth',2);
plot(xpcR,ypcR,'--k','LineWidth',2);
plot(xRT,yRT2,'-y','LineWidth',1);
plot(xRT1,yRT3,'-k','LineWidth',1);
end
    saveas(gcf, [ddir filesep subj filesep 'processed' filesep 'LGrid' filesep subj '_' tank], 'fig');
%     saveas(gcf, [ddir filesep subj filesep 'processed' filesep subj '_' tank 'b'], 'png');
    hold off
    figure('Name',[subj ' Electrode: ' tank ' Raw EEG']);
    x1 = 1:length(DAT);
    x1 = x1/(fs*60);
    plot(x1,DAT);
    if dobehave == 1 && gb == 1
    hold on
    plot(xpc,ypc,'--r','LineWidth',2);
    plot(xpcR,ypcR,'--k','LineWidth',2);
    plot(xRT,yRT2,'-y','LineWidth',1);
    plot(xRT1,yRT3,'-k','LineWidth',1);
    
    end
    saveas(gcf, [ddir filesep subj filesep 'processed' filesep 'LGrid' filesep subj '_' tank 'raw'], 'png');

end
end