function ram_freq_anesthesia
%
%   Compute power spectrum from data.
%
%


% parameters

ddir = '/Users/tnl/Desktop/C/data/eeg/'; % directory
subj = 'HUP119_i'; % subject
subjF = 'HUP119_i_24Jun16_0959'; % raw data name?
snipit = 6; % time, in minutes to clip at the beginning of the spectrogram
elec = 1; % Electrode number
tank = num2str(elec); % Electrode number
fs = 1000; % sampling rate (Hz) fter being resampled
nch = 1; % number of channels
etime = 23.5; % length of time in minutes
% cor = [3.92 15.233]; % correct range
% LWM = 10.519; % LWM (loss of memory) in minutes from events file (lfpoffset / samplerate / 60) samplerate from params.txt
% LOC = 15.756; % LOC in minutes from events file
% sr = 2.047999e+03;

sr = eegparams('samplerate',[ddir subj filesep 'lfp.noreref' filesep 'params.txt']); %obtain sample rate from params text file
timefix = sr*60; % conversion to minutes from time points?

events = [];
load([ddir subj filesep 'behavioral' filesep 'session_0' filesep 'events.mat']); % get behavioral data
fixEEGlfpfile_MAC % use this when working from the MAC
%fixEEGlfpfile; % changes the .lfpfile name to the correct string for windows
% eventStrategize loops through RECALL events nd checks if the resp eqials
% the target
[correctR incorrectR] = eventStrategize(events); % (correct range) find correct, incorrect, and pass (LOC) recall trials
cor = [correctR(1).lfpoffset/timefix correctR(size(correctR,2)).lfpoffset/timefix]; %offset of 1st and offset of last trial
LWM = incorrectR(1).lfpoffset/timefix; % Loss of working mem: end of 1st incorrect
LOC = correctR(size(correctR,2)).lfpoffset/timefix; % loss of consciousness: estmated based on final correct trial

% convert to minutes
for n=1:length(events)
    events(n).lfpoffset2 = events(n).lfpoffset/timefix;
end

% initialize counts?
ci = 1;
ii = 1;
rti = 1;
% sr = 2.047999e+03 %sample rate from

% quantify behavioral data and get data to plot
for i = 1:length(events)
    % if this is a string of 3 RECALL events,
    if i+2<length(events) && strcmp(events(i).type,'RECALL') && strcmp(events(i+2).type,'RECALL');
        % add "correct" field
        events(i).correct = strcmp(events(i).response,events(i).target);
        events(i+1).correct = strcmp(events(i+1).response,events(i+1).target);
        events(i+2).correct = strcmp(events(i+2).response,events(i+2).target);
        % add "percent correct"
        events(i).percentcorrect = (events(i).correct + events(i+2).correct + events(i+1).correct)/3*100;
        events(i+1).percentcorrect = events(i).percentcorrect;
        events(i+2).percentcorrect = events(i).percentcorrect;
        %         pc(ci) = events(i);
        %         pc(ci+1) = events(i+1);
        %         pc(ci+2) = events(i+2);
        % x points
        xpc(ci) = events(i).lfpoffset2;
        xpc(ci+1) = events(i+1).lfpoffset2;
        xpc(ci+2) = events(i+2).lfpoffset2;
        % y points
        ypc(ci) = events(i).percentcorrect;
        ypc(ci+1) = events(i+1).percentcorrect;
        ypc(ci+2) = events(i+2).percentcorrect;
        ci = ci + 3;
    end
    % if this is a string of 3 responses, calculate percecntages
    if i+4<length(events) && strcmp(events(i).type,'RESP') && strcmp(events(i+2).type,'RESP') && strcmp(events(i+4).type,'RESP');
        % was this correct
        events(i).correctRESP = strcmp(events(i).response,events(i).target);
        events(i+2).correctRESP = strcmp(events(i+2).response,events(i+2).target);
        events(i+4).correctRESP = strcmp(events(i+4).response,events(i+4).target);
        % percentage
        events(i).percentcorrectRESP = (events(i).correctRESP + events(i+2).correctRESP + events(i+4).correctRESP)/3*100;
        events(i+2).percentcorrectRESP = events(i).percentcorrectRESP;
        events(i+4).percentcorrectRESP = events(i).percentcorrectRESP;
        %         pcR(ii) = events(i);
        %         pcR(ii+1) = events(i+2);
        %         pcR(ii+2) = events(i+4);
        % x and y points
        xpcR(ii) = events(i).lfpoffset2;
        xpcR(ii+1) = events(i+2).lfpoffset2;
        xpcR(ii+2) = events(i+4).lfpoffset2;
        ypcR(ii) = events(i).percentcorrectRESP;
        ypcR(ii+1) = events(i+2).percentcorrectRESP;
        ypcR(ii+2) = events(i+4).percentcorrectRESP;
        ii = ii + 3;
    end
    
    % if this is a single response (in string of 3), calculate RT
    if i+1<length(events) && strcmp(events(i).type,'SOUND') && strcmp(events(i+1).type,'RESP')
   % rt is difference between the response (i+1) and the end of the sound
        events(i).RT = events(i+1).lfpoffset-events(i).lfpoffset;
        % x and y points for plotting
        xRT(rti) = events(i).lfpoffset2;
        yRT(rti) = events(i).RT/sr;
        rti = rti + 1;
    end
    
end

yRT2 = yRT/max(yRT)*100; % percent maximal RT. max RT is when the subject doesnt respond

events = events(1); % take only the first event to look at entire time sequence.
% events.lfpfile = [ddir subj filesep 'lfp.noreref' filesep subjF];


% data
freq = logspace(log10(0.5),log10(200),100)'; % frequencies (Hz)
winsz = 5; % window size for spectral estimates (s) - determines frequency resolution
% try load([ddir subj filesep 'processed' filesep subj '_' tank '.mat'],'-mat');
% catch
%     DAT = [];
% this will return raw LFP traces, also downsamples and filters
DAT = an_getlfp_ms_wrapper(elec,events,etime*60*1000,0,0,[0 200],'bandpass',1);

%     DAT = an_getlfp_ms_wrapper(elec,events,etime*60*1000,0,0);


disp('loading data...');

lastfullcycle = floor(length(DAT)/nch)*nch;
DAT(lastfullcycle+1:end) = [];
DAT = reshape(DAT,nch,[])';
N = size(DAT,1); % samples
S = fix(N/fs)-(winsz-1); % samples to seconds, - the window size. total seconds to plot?

% spectra
PS = zeros(nch,S,length(freq)); %power spectrum, channels by seconds by freq
fprintf('computing spectra...');
for ii = 1:S % compute in 'winsz'-sec windows in 1-sec steps
    fprintf('%02d%%',fix(ii/S*100));
    ind = round((ii-1)*fs+1:(ii+winsz-1)*fs); % get indicies of window
    for jj = 1:nch
        PS(jj,ii,:) = periodogram(DAT(ind,jj),hamming(length(ind)),freq,fs); % calculated PSD % TODO: remove mean before spectral analysis?
    end
    if ii~=S, fprintf('\b\b\b'); end
end
fprintf('\n');

% smoothing
M = fix(S/60); % number of minutes?
PM = zeros(nch,size(PS,3),M);
disp('smoothing spectra...');
for ii = 1:M
    ind = (ii-1)*60+1:ii*60; % indicied of  minutes
    for jj = 1:nch
        PM(jj,:,ii) = trimmean(squeeze(PS(jj,ind,:)),10); % robust mean for every 60-sec (removes upper and lower 'winsz' seconds of each minute of data)
    end
end
x = 1:M;



% save
%     cd_mkdir([ddir subj filesep 'processed' filesep tank]);
%save([ddir subj filesep 'processed' filesep subj '_' tank],'PM','x','freq','DAT','PS','sr','xpc','ypc','xpcR','ypcR','xRT','yRT','yRT2');
% end

% pre-plot
ind = find(x>snipit); % clipping beginning of spectrogram
x = x(ind);
PM = PM(:,:,ind);
for ich = 1:size(PM,1)
    PM(ich,:,:) = smooth2(squeeze(PM(ich,:,:)),2,1); % (data, dx, dy)
    PM(ich,:,:) = zscore(squeeze(PM(ich,:,:))')'; % zscore
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
    set(gca,'CLim',[-3 4]);
    hc = colorbar; set(get(hc,'YLabel'),'String','power');
    %     line([LWM LWM],[0 200],'LineWidth',4,'Color','k')
    %     line([LOC LOC],[0 200],'LineWidth',4,'Color','r')
    %     line(cor,[1 1],'LineWidth',2,'Color','r')
    hold on
    plot(xpc,ypc,'--r','LineWidth',2); % correct recal
    plot(xpcR,ypcR,'--k','LineWidth',2); % correct resp
    plot(xRT,yRT2,'-y','LineWidth',1); % RT% mx
    %saveas(gcf, [ddir subj filesep 'processed' filesep subj '_' tank], 'fig');
end