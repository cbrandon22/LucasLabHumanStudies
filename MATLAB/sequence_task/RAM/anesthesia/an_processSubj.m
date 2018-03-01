function[] = an_processSubj(subj,sess)
%This function is a wrapper to process a single subject's data. 
%1) It makes events using an_makeEvents.m
%2) it processes raw neuralynx files using an_downsampleAndSaveNlxData
%3) it aligns the behavioral clock with the nlx clock using
%an_alignNeuralynxFile and updates the events file appropriately (adds the
%following fields: lfpfile, lfpoffset (samples), timesfile, timesOffset
%(ms))

%Inputs
%subj.......'NR_062415'
%sess.......'post'    
%lfpFile.....(e.g., 'TJ061_30Apr13_1240'; lfp file name to locate sync
             %pulses). Output from downsampleAndSaveNlxdata
%% get dirs etc
[subjInfo,dataDir] = an_subjInfo(subj);
sessInd = find(strcmp(subjInfo.sess,sess));
behDir = fullfile(dataDir,subj,'behavioral',['session_' sess]);
%% process behavioral data (Create events structure)
an_makeEvents(subj,sess);
an_fixEEGLog('eeg.eeglog','eeg.eeglog.up');

%% process neural data (save downsampled nlx files)
lfpFile = an_downsampleAndSaveNlxData(subj,subjInfo.rawDataPath{sessInd},subjInfo.logpeg);

%% filenames
lfpFullFile = fullfile(dataDir,subj,'lfp.noreref',lfpFile);
unitFullFile = fullfile(dataDir,subj,'units',subjInfo.rawDataPath{sessInd},'times_CSC'); %this is a place holder for future use

%% list of events to align
listOfEvents{1}=fullfile(behDir,'events.mat');

%% run align
an_alignNeuralynxFile(lfpFullFile,behDir,listOfEvents,unitFullFile)
