function ccdtCorrPrep
% function ccdtCorrPrep
%   Neural correlates of CCDT: data preparation. Save spectrogram aligned
%   on color change event for each rewarded trials. Also save reaction time
%   and delay time information and error trial information.
%
%   Valid for sessions starting with 20151009 (9 calibration targets & eye
%   signals converted to pixels).
%
%   DR 03/2016

% parameters
ddir = '/Users/tnl/Desktop/'; % data directory
subj = 'fidel'; % subject
dtank = 'fidel_VPLT_DT1_110415'; % tank
dblock = 1; % block
ch = 1:16; % channels
twin = [-2000 1500]; % time window relative to color change (ms)
twine = [-500 2500]; % time window realtive to cue for error trials (ms)
dns = 50; % downsample factor
bpf = [1 200]; % bandpass filter (Hz)
freq = logspace(log10(1),log10(200),50)'; % frequencies (Hz)
wavl = 'cmor1-1'; % complex wavelet
pdir = [ddir subj filesep 'processed' filesep 'ccdtCorr' filesep]; % prepared data directory

% constants
fs = 24414.06; % sampling rate (Hz)
scales = centfrq(wavl)./(freq*(1/(fs/dns)));
Nf = length(freq); % number of frequencies

% behavioral data
cd([ddir subj filesep dtank filesep 'Block-' num2str(dblock)]);
fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(19) '.sev'],'r');
fseek(fid,40,'bof');
CMD = fread(fid,inf,'single')/1e6; % plotting command
fclose(fid);
fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(20) '.sev'],'r');
fseek(fid,40,'bof');
PD = fread(fid,inf,'single')/1e6; % photodetector signal
fclose(fid);
iccc = find(CMD(1:end-1)>=211 & CMD(1:end-1)<=219 & CMD(2:end)==220); % color change command
irp = find(PD(1:end-1)<0.8 & PD(2:end)>0.8); % respond plot
icc = zeros(size(iccc));
for ii = 1:length(iccc)
    [~,ind] = min(abs(irp-iccc(ii)));
    icc(ii) = irp(ind); % color change (actual)
end
t = (1:length(CMD))/fs;
talign = t(icc); % time of color change (s)
targ = CMD(iccc-100); % target code (211-219)
iin = find(CMD(1:end-1)==220 & (CMD(2:end)==210 | CMD(2:end)==0)); % CCDT response (or timeout)
try
    RT = (iin-icc)/fs*1000; % reaction time
catch
    disp('something wrong with behavioral data'); return;
end
ifp = find(PD(1:end-1)<0.2 & PD(2:end)>0.2); % fixation plot
iff = zeros(size(icc));
for ii = 1:length(icc)
    [~,ind] = min(abs(ifp-icc(ii)));
    iff(ii) = ifp(ind); % color change (actual)
end
DT = (icc-iff)/fs*1000; % delay time
ibad = find(RT<=0 | RT>=1000 | DT<0); % only include rewarded trials
talign(ibad) = []; RT(ibad) = []; DT(ibad) = []; targ(ibad) = [];
Nt = length(talign); % number of trials
if isempty(RT)
    disp('no behavioral data'); return;
end
ipr = find(CMD(1:end-1)>=211 & CMD(1:end-1)<=219 & (CMD(2:end)==210 | CMD(2:end)==0)); % response before color change (assumes 9 calibration targets)
ifx = zeros(size(ipr));
for ii = 1:length(ipr)
    [~,ind] = min(abs(ifp-ipr(ii)));
    ifx(ii) = ifp(ind);
end
DTe = (ipr-ifx)/fs*1000; % delay time on error trials
taligne = t(ifx); % time of cue on error trials

% gaze data
fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(17) '.sev'],'r');
fseek(fid,40,'bof');
eye0 = fread(fid,inf,'single')/1e6; % eye X (pixels) divide by 1 million for Volts (system records in microvolts)
fclose(fid);
fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(18) '.sev'],'r');
fseek(fid,40,'bof');
eye1 = fread(fid,inf,'single')/1e6; % eye Y (pixels)
fclose(fid);
Sindex = saccadeFinder(eye0,eye1,fs);
eyeX = eye0;
eyeY = eye1;
ind = 1:Sindex(1,1);
eyeX(ind) = median(eye0(ind));
eyeY(ind) = median(eye1(ind));
for ii = 2:size(Sindex,1) % remove artifacts/blinks by computing robust gaze position between saccades (assumes saccadeFinder is robust to artifacts)
    ind = Sindex(ii-1,2):Sindex(ii,1); 
    eyeX(ind) = median(eye0(ind));
    eyeY(ind) = median(eye1(ind));
end
ind = Sindex(end,2):length(eyeX);
eyeX(ind) = median(eye0(ind));
eyeY(ind) = median(eye1(ind));
eyeX = downsample(eyeX,dns);
eyeY = downsample(eyeY,dns);

% neural data
cd([ddir subj filesep dtank filesep 'Block-' num2str(dblock)]);
[blp,alp] = butter(2,bpf/(fs/2)); 
for ich = 1:length(ch)
    fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(ch(ich)) '.sev'],'r');
    fseek(fid,40,'bof');
    dat = fread(fid,inf,'single');
    fclose(fid);
    dat = filtfilt(blp,alp,dat); % anti-aliasing filter
    dat = downsample(dat,dns); % downsample
    if mean(dat)==0
        disp('no neural data'); continue;
    end
    
    % wavelet transform
    Pxx = cwt(dat,scales,wavl); % complex wavelet coefficients
    
    % window
    if ich==1
        fs = fs/dns;
        t = (1:length(dat))/fs;
        salign = zeros(Nt,1);
        for ii = 1:Nt
            [~,salign(ii)] = min(abs(t-talign(ii)));
        end
        swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
        indwin = salign*ones(1,length(swin)) + ones(Nt,1)*swin;
        ibad = find(indwin(:,1)<1 | indwin(:,end)>length(dat));
        indwin(ibad,:) = []; RT(ibad) = []; DT(ibad) = []; targ(ibad) = []; Nt = size(indwin,1);
        eyeXcc = eyeX(indwin);
        eyeYcc = eyeY(indwin);
        saligne = zeros(length(taligne),1);
        for ii = 1:length(taligne)
            [~,saligne(ii)] = min(abs(t-taligne(ii)));
        end
        swine = round(twine(1)/1000*fs):round(twine(2)/1000*fs);
        indwine = saligne*ones(1,length(swine)) + ones(length(taligne),1)*swine;
        ibad = find(indwine(:,1)<1 | indwine(:,end)>length(dat));
        indwine(ibad,:) = []; DTe(ibad) = [];
        t = swin/fs*1000;
        te = swine/fs*1000;
    end
    Pcc = zeros(Nt,Nf,size(indwin,2));
    Pcce = zeros(size(indwine,1),Nf,size(indwine,2));
    for ii = 1:Nf
        cP = Pxx(ii,:)';
        Pcc(:,ii,:) = cP(indwin);
        Pcce(:,ii,:) = cP(indwine);
    end 
    Xcc = dat(indwin);
    Xcce = dat(indwine);
    
    % save
    save([pdir dtank '_' num2str(dblock) '_' num2str(ch(ich)) '.mat'],'RT','DT','DTe','targ','Pcc','Pcce','Xcc','Xcce','eyeXcc','eyeYcc','t','te','freq','-mat');
end
