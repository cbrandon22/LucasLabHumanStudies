% This script was written to test the effect of different
% Neuralynx settings on the power spectral density. ddir is the path to the
% CSC files from each sample setting. CB -6/13/2018

% First set of 4 recordings
% ddir={'/Volumes/Lucas Drive/exports/nlxTest/1kHFF_11khz/2018-06-13_13-01-38';...
% 	'/Volumes/Lucas Drive/exports/nlxTest/1kHFF_32khz/2018-06-13_15-47-09';...
% 	'/Volumes/Lucas Drive/exports/nlxTest/9kHFF_11khz/2018-06-13_16-01-37';...
% 	'/Volumes/Lucas Drive/exports/nlxTest/9kHFF_32khz/2018-06-13_14-08-33';};

% Second set of 4 recordings
% ddir={'/Volumes/Lucas Drive/exports/nlxTest/1kHFF_11khz/2018-06-13_16-27-29';...
% 	'/Volumes/Lucas Drive/exports/nlxTest/1kHFF_32khz/2018-06-13_16-56-59';...
% 	'/Volumes/Lucas Drive/exports/nlxTest/9kHFF_11khz/2018-06-13_17-10-58';...
% 	'/Volumes/Lucas Drive/exports/nlxTest/9kHFF_32khz/2018-06-13_16-42-15';};

% Third set of 4 recordings (using faraday cage)
% ddir={'/Volumes/Lucas Drive/exports/nlxTest/1kHFF_11khz/2018-06-14_08-55-43';...
% 	'/Volumes/Lucas Drive/exports/nlxTest/1kHFF_32khz/2018-06-14_09-03-25';...
% 	'/Volumes/Lucas Drive/exports/nlxTest/9kHFF_11khz/2018-06-14_12-46-35';...
% 	'/Volumes/Lucas Drive/exports/nlxTest/9kHFF_32khz/2018-06-14_12-54-31';};

% Forth set of 2 recordings (EMU vs OR - real patient)
ddir={'/Volumes/HumanStudies/HumanStudies/oddball/eeg/HUP171_emu/raw/2018-06-25_16-52-58';...
    '/Volumes/HumanStudies/HumanStudies/oddball/eeg/HUP171_emu/raw/1kHFF';...
	'/Volumes/HumanStudies/HumanStudies/oddball/eeg/HUP171_i/raw/2018-06-26_07-19-07';...
    '/Volumes/Lucas Drive/exports/nlxTest/9kHFF_32khz/2018-06-14_12-54-31'};

freq = [0.5 10000];
lowpassfreq = 1000;
subjList = {'HUP155_i'};
subj_fn='/Volumes/HumanStudies/HumanStudies/oddball/eeg/HUP155_i/raw/2017-12-15_19-39-14/CSC1.ncs';
[data,info,~,~] = load_ncs2(subj_fn);
data = double(data);
fs = info.actualSampleRate;
[pxx_cat,f_cat] = pwelch(data,5*fix(fs),1*fix(fs),logspace(log10(freq(1)),log10(freq(2)),750),fs); % 5-sec windows, 1-sec overlap, 1000 logorithmically-spaced power estimates
pxx_lp_cat = [];
f_lp_cat = [];
for i=1:length(ddir)
    elec = 72;                          % electrode
    subj = strsplit(ddir{i},'/');
    subj = subj{end-2};
    if strcmp(subj,'nlxTest'),elec=4;end
    if i==2,subj = [subj '_1kHFF'];end
    subjList = [subjList,subj];
    % load data
    cd(ddir{i});
    [data,info,~,~] = load_ncs2(['CSC' num2str(elec) '.ncs']);
    %[data,info,~,~] = load_ncs2(['MW' num2str(elec) '.ncs']);
    data = double(data);
    fs = info.actualSampleRate;
    data_lp = buttfilt(data,lowpassfreq,fs,'low',4);
    
    % Welch's PSD
    [pxx,f] = pwelch(data,5*fix(fs),1*fix(fs),logspace(log10(freq(1)),log10(freq(2)),750),fs); % 5-sec windows, 1-sec overlap, 1000 logorithmically-spaced power estimates
    pxx_cat = [pxx_cat;pxx];
    f_cat = [f_cat;f];
    
    % Welch's PSD
    [pxx_lp,f_lp] = pwelch(data_lp,5*fix(fs),1*fix(fs),logspace(log10(freq(1)),log10(freq(2)),750),fs); % 5-sec windows, 1-sec overlap, 1000 logorithmically-spaced power estimates
    pxx_lp_cat = [pxx_lp_cat;pxx_lp];
    f_lp_cat = [f_lp_cat;f_lp];
end
subjList(end) = {'OR Saline'};
% plot
figure('Name','HUP171 Test','NumberTitle','off','Units','normalized','Position',[1/4 1/4 1/2 1/2],'Color','w');
hold on;
for i=2:3%size(pxx_cat,1)
    plot(f_cat(i,:),pxx_cat(i,:));
end
axis tight; set(gca,'Box','off','XScale','log','YScale','log');
xlabel('frequency (Hz)'); ylabel('power (\muV^2/Hz)'); title('5-sec windows, 1-sec overlap, 750 log-spaced power estimates HUP155/HUP171/saline recordings','FontSize',20);
legend(subjList(2:3),'Interpreter', 'none','FontSize',20);
keyboard
% Low pass filter, downsample, recalculate spectrum
figure('Name','HUP171 low-pass filter','NumberTitle','off','Units','normalized','Position',[1/4 1/4 1/2 1/2],'Color','w');
hold on;
for i=1:size(pxx_lp_cat,1)
    plot(f_lp_cat(i,:),pxx_lp_cat(i,:));
end
axis tight; set(gca,'Box','off','XScale','log','YScale','log');
xlabel('frequency (Hz)'); ylabel('power (\muV^2/Hz)'); title('5-sec windows, 1-sec overlap, 750 log-spaced power estimates HUP171/saline recordings','FontSize',20);
legend(subjList(2:end),'Interpreter', 'none','FontSize',20);