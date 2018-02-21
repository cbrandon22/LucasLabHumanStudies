clear

subj           = 'TJ061';
pathToBehavDir = 'prob_sel2/TJ061_Olivia/session_0_train/';
fileExt        = 'TJ061_30Apr13_1240';
timesfile = '/data/eeg/TJ061/raw/2013-04-30/nlx/times_CSC';
evCellArray{1} = 'events.mat'; 
%evCellArray{2} = 'MATH_events.mat';
%evCellArray{3} = 'RECOG_events.mat';

%-----------------------------------------
%-----------------------------------------
dataDir     = fullfile('/data/eeg/',subj);
lfpDir      = fullfile(dataDir,'lfp.noreref');
behDir      = fullfile(dataDir,'behavioral');
fullFileExt = fullfile(lfpDir,fileExt);
fullBehPath = fullfile(behDir,pathToBehavDir);

%-------------------------------------------
for k=1:length(evCellArray)
  evCellArray{k}=fullfile(fullBehPath,evCellArray{k});
end

%-------------------------------------------
alignNeuralynxFile(fullFileExt,fullBehPath,evCellArray,timesfile)