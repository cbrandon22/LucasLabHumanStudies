subj = 'HUP155_i';
ddir = '/Volumes/HumanStudies/HumanStudies/oddball/eeg'; %path to folder containing subjects
cscdir = fullfile(ddir,subj,'raw','2017-12-15_19-39-14');
win_ms = 300000;
load(fullfile(ddir,subj,'processed/sessInfo.mat'));
load(fullfile(ddir,subj,'processed/parameters.mat'));
if exist(fullfile(ddir,subj,'processed/manual_bad_channels.mat'),'file')==2
    load(fullfile(ddir,subj,'processed/manual_bad_channels.mat'));
end

win = round(win_ms*(srate/1000));
startWin = srate*120;

good_elecInfo = elecInfo(~ismember(cell2mat(elecInfo(:,1)),bad_channels),:);

data = ones(size(good_elecInfo,1),win)*NaN;
for chan = 1:size(good_elecInfo,1)
    dat = look(events(1).lfpfile,good_elecInfo{chan,1},[],1)';
    data(chan,:) = dat(startWin:startWin+win-1);
    %[dat,info,~,~] = load_ncs2([cscdir '/' good_elecInfo{chan,2} '.ncs']);
    %dat = double(dat);
end
keyboard

figure
hold on
plot_ms = 100;
t = linspace(0,plot_ms,plot_ms*(srate/1000));
for chan = 1:size(good_elecInfo,1)
    plot(t,data(chan,1:round(plot_ms*(srate/1000))));
    ylabel('Voltage');
    xlabel('Time (ms)');
end

freq = [0.1 10000];
nw = 3;
pxx = pmtm(data(1,:),3,linspace(freq(1),freq(2),20000),srate);
plot(pxx)
set(gca,'XScale','log','YScale','log');

figure;
freq = [0.1 10000];
nw = 3;
pxx = pmtm(data(1,:),3,logspace(log10(freq(1)),log10(freq(2)),10000),srate);
plot(pxx)
set(gca,'XScale','log','YScale','log');

