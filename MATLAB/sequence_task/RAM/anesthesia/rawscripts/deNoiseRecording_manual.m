function [] = deNoiseRecording_manual(rawDir,rawFile,sr,fmin,fmax)
% This function loads the data, filters it as wave_clus does, and plots the
% entire session, and then keyboards. The user can then entire the time
% periods to eliminate. It then saves out a data file with the same path
% but an '_denoised_man' appended to it.
% Inputs:
% rawDir            where the raw data is saved (eg '/data5/eeg/CABG061/units')
% rawFile           name of the recorded file (eg 'CABG061-LSN-210807LSTNTrack1-1')
% sr                samplerate of the recording
% fmin (optional)   default minimum frequency
% fmax (optional)   default maximum frequency

% To Do:
% ( ) make it compatible with ncs recordings
% ( ) write a gui for this
% written by ashwin ramayya 12-11-2013

%%% REMOVE
sr = 25000;

% Load data
load(fullfile(rawDir,rawFile))
x = data;

%filter params
if ~exist('fmin','var') || isempty(fmin)
    fmin = 400;
end
if ~exist('fmax','var') || isempty(fmax)
    fmax = 3000;
end

% filter
[b,a] = ellip(2,0.1,40,[fmin fmax]*2/(sr));
data_filt =  filtfilt(b,a,x);

h(1) = swagFig([-0.0020    0.7959    0.8000    0.1553]);
plot(data_filt);
swagAxes(gca,12)
cd(rawDir)
keyboard

%%%%% START HERE

% denoise via threshold
deNoise_lims = [-200 200];
data_filt_deNoised = deNoise_local(data_filt,deNoise_lims);
h(1) = swagFig([-0.0020    0.7959    0.8000    0.1553]);
plot(data_filt_deNoised);
swagAxes(gca,12)


% set a manual thresholdl % plot 500 ms clips
manThresh_std_sua = 4; manThresh_sua = manThresh_std_sua*nanstd(data_filt_deNoised);
plotRandClips_local(data_filt_deNoised,sr,manThresh_sua)

manThresh_std_mua = 2.5; manThresh_mua = manThresh_std_mua*nanstd(data_filt_deNoised);
plotRandClips_local(data_filt_deNoised,sr,manThresh_mua)

% save
cd(rawDir)
% re-write data vector
data(isnan(data_filt_deNoised)) = 0;
save([rawFile '_denoised_manual.mat'],'data','manThresh_mua',...
    'manThresh_std_mua','manThresh_sua','manThresh_std_sua','deNoise_lims','fmin','fmax','sr')


function data_filt_deNoised = deNoise_local(data_filt,data_lims)
samp2rm = find(data_filt<data_lims(1) | data_filt>data_lims(2));
samples_before = 10;
samples_after = 10;
samp2rm = samp2rm((samp2rm > samples_before) & (samp2rm < length(data_filt) - samples_after));
data_filt_deNoised = data_filt;
for t = -samples_before:samples_after %tricky!!
    data_filt_deNoised(samp2rm+t) = nan;
end
%data_filt_deNoised(samp2rm) = nan;

function plotRandClips_local(data,sr,thresh, thresh2)
if ~exist('thresh2','var') || isempty(thresh2)
    thresh2 = thresh;
end
swagFig([0.7504    0.0082    0.2234    0.9430])
dur_samp = (500*sr)/1000;
for i = 1:10
    subplot(10,1,i)
    rand_start = randi([dur_samp length(data)-(dur_samp+1)],[1 1]); 
    plot(data(rand_start:rand_start+dur_samp));swagAxes(gca,10)
    hold all
    plot(get(gca,'xlim'),[-thresh -thresh],'r')
    plot(get(gca,'xlim'),[thresh thresh],'r')
    plot(get(gca,'xlim'),[-thresh2 -thresh2],'--r')
    plot(get(gca,'xlim'),[thresh2 thresh2],'--r')
end