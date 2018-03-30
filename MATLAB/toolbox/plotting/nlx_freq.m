function nlx_freq
% function nlx_freq
% 	Power spectral density of one raw signal recorded on Neuralynx system.
%
%   DR 10/2017

% parameters
dirs = le_dirs('oddball');
ddir = fullfile(dirs.data,'eeg/');
%ddir = 'E:\TNL_Data\oddball\eeg\';  % data directory
subj = 'HUP149_e';                  % subject
elec = 16;                          % electrode
freq = [0.5 500];                   % frequency range (Hz)

% load data
cd([ddir subj '/raw']);
fl = dir('*-*');
cd(fl(1).name);
[data,info,~,~] = load_ncs2(['CSC' num2str(elec) '.ncs']);
data = double(data);
fs = info.actualSampleRate;

% Welch's PSD
[pxx,f] = pwelch(data,10*fix(fs),1*fix(fs),logspace(log10(freq(1)),log10(freq(2)),100),fs); % 10-sec windows, 1-sec overlap, 100 logorithmically-spaced power estimates

% plot
figure('Name',subj,'NumberTitle','off','Units','normalized','Position',[1/4 1/4 1/2 1/2],'Color','w');
plot(f,pxx);
axis tight; set(gca,'Box','off','XScale','log','YScale','log');
xlabel('frequency (Hz)'); ylabel('power (\muV^2/Hz)'); title([subj ' elec ' num2str(elec)]);
