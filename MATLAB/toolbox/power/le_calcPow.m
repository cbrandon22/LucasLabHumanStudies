function[subjPow,config] = le_calcPow(subj,sessList,elecNum,task,evLabel,configNum)
%% function[] = le_calcPow(configNum,elecNum,subj,sess,task,evLabel)
%This function loads power for a particular electrode or electrode pair
%associated with a particular subject's events from a particular task.
%Filters by session so that power from relevant sessions are returned. 
%If power is being computed for the first time, it saves it in the scratch
%directory defined in le_dirs

%Inputs
%subj  = 'UP034'; 
%sessList = {'0'}; % (optional; if empty, it loops through all sessions seperately all events); 
                %Note: must be a cell array of strings
%elecNum = [1 2] (for bipolar contacts)
%task = 'pyFR'; %the label of the task dir in /data/events
%evLabel = 'math' (optional: default is 'events', but allows you to load
          %math or recog from pyFR if needed)

%Outputs
%subjPow        concatonated power for each events (freq x time x event); 
                %see config file to get time and freq windows
%Written by AGR 09-23-2013
%% set dirs and parse inputs
dirs = le_dirs(task);
if ~exist('evLabel','var') || isempty(evLabel)
    evLabel = 'events';
end
if ~exist('configNum','var') || isempty(configNum)
    configNum = 1;
end

%% load config
config = le_config_calcPow(configNum,task);

%% go to pow directory
powDir = fullfile(dirs.scratch,'POWER',num2str(config.configID),subj,[task '_' evLabel]);
cd_mkdir(powDir);
subjPow = [];

%% load events
%load events
[events] = le_loadEvents(subj,sessList,task,evLabel);
%unique([events.session])
%% if ~exist sessList, create it here
if ~exist('sessList','var') || isempty(sessList)
    sessList = unique({events.session});
end

%% loop through sessList
for i = 1:length(sessList)
    sess = sessList{i};
    cd_mkdir(fullfile(powDir,sess))
    if config.bipolFlag
        powFile = [num2str(elecNum(1)) '-' num2str(elecNum(2)) '.mat'];
    else
        powFile = [num2str(elecNum) '.mat'];
    end
    
    %try and load the power file
    if exist(powFile,'file')
        load(powFile)
    else 
        disp([num2str(elecNum) '-' sess]);
        %filter events
        sessEv = events(strcmp({events.session},sess));
        %calculate power for this session and save
        pow = single(calcPow(elecNum,sessEv,config));
        save(powFile,'pow')
    end
    
    %concatonate pow
    subjPow = cat(3,subjPow,pow);
end

%% functions below
function[pow] = calcPow(elecNum,events,config)
%this function calculates power for a particular session
%gete_ms_wrapper
if ~config.nlxFlag
    [EEG] = gete_ms_wrapper(elecNum,events, ...
               config.durationMS,config.offsetMS,config.bufferMS,...
               config.filtfreq, config.filttype, config.filtorder);
else
    [EEG] = an_getlfp_ms_wrapper(elecNum,events, ...
               config.durationMS,config.offsetMS,config.bufferMS,...
               config.filtfreq, config.filttype, config.filtorder);
end
%eeg2pow 
[pow] = eeg2PowAndPhase(EEG,config.freQ,config.bufferMS,config.logBefMeanFlag,[],config.waveNum); 
[pow] = pow2timeFreqBins(pow,config.timeWin,config.timeStep,config.freQ,config.freqBins,...
    config.logAfterMeanFlag);  
if config.zFlag==1
    % get baseline
    [baseEv]=get_baseline_random(events,config.baseSecBetEv,config.baseJitt);
    if ~config.nlxFlag
        [baseEEG] = gete_ms_wrapper(elecNum,baseEv, ...
           config.basePriorMS+config.basePostMS,...
           -config.basePriorMS,config.bufferMS);
    else
        [baseEEG] = an_gete_ms(elecNum,events, ...
               config.durationMS,config.offsetMS,config.bufferMS,...
               config.filtfreq, config.filttype, config.filtorder);
    end
    [rawBasePow] = eeg2PowAndPhase(baseEEG,config.freQ,config.bufferMS,config.logBefMeanFlag,[],config.waveNum); 
    meanRawBasePow = mean(mean(rawBasePow,2),3);
    stdRawBasePow = std(mean(rawBasePow,2),[],3);
    [basePow] = pow2timeFreqBins(rawBasePow,config.timeWin,config.timeStep,config.freQ,config.freqBins,...
               config.logAfterMeanFlag);  
    basePow = shiftdim(nanmean(basePow,2),2)';
    meanBasePow = nanmean(basePow,2);
    stdBasePow = nanstd(basePow,[],2);

    %z score
    pow = zScorePow(pow,meanBasePow,stdBasePow);
    pow = single(pow);
end


function[zPowMat] = zScorePow (powMat,meanBasePow, stdBasePow)
%This function z-scores a power matrix by a baseline power dist, frequency, by frequency
%Inputs 
%powMat                         %freq x time x trials
%meanBasePow                    %mean pow for each freq
%stdBasePow                     %std of pow at each freq
%Output
%zPowMat                        %zScored power matrix

%Written by AGR 11-04-2012

zPowMat = (powMat - repmat(meanBasePow,[1 size(powMat,2), size(powMat,3)]))...
    ./repmat(stdBasePow,[1 size(powMat,2), size(powMat,3)]);

function[eeg,resampfreq] = gete_ms_wrapper(elecNum,events, ...
               durationMS,offsetMS,bufferMS,filtfreq,filttype,filtorder)
%This wrapper returns raw eeg traces for a set of events. 
%Inputs
%elecNum        %electrode number to load. Can be a single channel (eg,[1]) or a set
                %of channels (for bipolar, eg,[1 2]). see getBipolarSubjElecs in
                %eeg_toolbox
%events         events structure (with 'eegfile', and 'eegoffset' fields).
                %see ps2_makeEvents and alignTool
%optional:
%durationMS     %how long to clip
%offsetMS       %where to start relative to mstime
%bufferMS       %buffer to add to beginning and end of clips
%filtfreq       %filters out particular frequencies (e.g., [59.5 60.5]);
                %see buttfilt
%filtorder      %e.g, 4 (see buttfult)

%Note:
%resampleFreq (last arg of gete) is set to 1000 so that gete returns eeg in ms.
                
%Written by AGR 10-31-12

%Parse Inputs
%
if ~exist('durationMS','var') || isempty(durationMS)
    durationMS = 1500;
end
if ~exist('offsetMS','var') || isempty(offsetMS)
    offsetMS = -500;
end    
if ~exist('bufferMS','var') || isempty(bufferMS)
    bufferMS = 0;
end
if ~exist('filtfreq','var') || isempty(filtfreq)
    filtfreq = [];
end 
if ~exist('filttype','var') || isempty(filttype)
    filttype = [];
end
if ~exist('filtorder','var') || isempty(filtorder)
    filtorder = [];
end
% if ~exist('resampleFreq','var') || isempty(resampleFreq)
%     resampleFreq = [];
% end

% Update durationMS and offsetMS with buffer
%(Check this: buffer input to gete_ms doesnt seem to work when not
%using a filter)
durationMS = durationMS + 2*bufferMS;
offsetMS = offsetMS-bufferMS;

%check whether it is a bipolar montage or not
if size(elecNum,2) == 1 % not bipolar, proceed normally
   [eeg,resampfreq] = gete_ms(elecNum, events,durationMS,offsetMS,bufferMS,...
            filtfreq,filttype,filtorder,1000);
    
elseif size(elecNum,2) == 2 % then it is bipolar
    %First, rename events.eegfile to point to ,noreref
%     for k=1:length(events)
%          events(k).eegfile=regexprep(events(k).eegfile,'eeg.reref','eeg.noreref');
%     end
    
    %Second, gete from both channels and subtract them from each other 
    [eeg1,resampfreq] = gete_ms(elecNum(1), events,durationMS,offsetMS,bufferMS,...
        filtfreq,filttype,filtorder,1000);
    [eeg2,resampfreq] = gete_ms(elecNum(2), events,durationMS,offsetMS,bufferMS,...
        filtfreq,filttype,filtorder,1000);
    
    eeg = eeg1 - eeg2;
end
function[powEEG,phaseEEG] = eeg2PowAndPhase(eeg,freQ,buffToRem,logBefMeanFlag,...
sr,waveNum)
%This function calculates phase and power values
%for a given eeg matrix. This is a wrapper for multiphasevec2
%Inputs 
%eeg                    %Matrix of eeg values. Rows = Trials
                        %columns = samples
%freQ                   %freQ of interest (def: 2-95 Hz, log spaced)
%buffToRem              %buffer to remove at the beginning and end 
%logBefMeanFlag         %If this is set to 1, then the function returns
                        %log10(power)
%sr                     %sample rate (def: 1000, to work with gete_ms)
%waveNum                %wave number for multiPhaseVec (default = 5)

%Written by AGR 11-1-12

%Identify frequencies
if ~exist('freQ','var') || isempty(freQ)
    freQ=logspace(log10(2), log10(95),30);
end
if ~exist('buffToRem','var') || isempty(buffToRem)
    buffToRem=[];
end
if ~exist('logBefMeanFlag','var') || isempty(logBefMeanFlag)
    logBefMeanFlag=1;
end
%Identify samp rate
if ~exist('sr','var') || isempty(sr)
    sr=1000;
end
%set wave number
if ~exist('waveNum','var') || isempty(waveNum)
    waveNum=7;
end

%Initialize variables
powEEG = single(nan(length(freQ),size(eeg,2),size(eeg,1)));
phaseEEG = powEEG;

%Loop through  trials
for i = 1:size(eeg,1)
    [phaseEEG(:,:,i), powEEG(:,:,i)]=multiphasevec2(freQ, eeg(i,:), sr,waveNum);
end

%Remove buffer at beginning and at the end
if ~isempty(buffToRem)
    buffInd = [1:buffToRem, (size(eeg,2)-buffToRem+1):size(eeg,2)];
    powEEG (:,buffInd,:) = [];
    phaseEEG (:,buffInd,:) = [];
end

% Log if logBefMean Flag is 1
if logBefMeanFlag==1
    powEEG = log10(powEEG);
end
function [binPowMat] = pow2timeFreqBins(powMat,timeWin,timeStep,freQ,freqBins,logAfterMeanFlag)     
%This function calculates average power within frequency and time bins.
%Inputs
%powMat                 %matrix of power values
%timeWin                %length of sliding window
%timeStep               %ms by which to advance the sliding window
%freQ                   %frequencies to use for wavelets
%freqBin                %ends of each freq bin (e.g. [ 40 95] means
                        %two bins, 2 - 40 and 40 - 95
%logAfterMeanFlag       %if this is set to 1, it calculates
                        %a log of binPowMat

%Outputs
%binPowMat              %power binned in frequency and time
%

%Written by AGR 11-3-12

% create time bins
timeBins = timeWin:timeStep:size(powMat,2);

%Initialize power mat 
binPowMat = nan([length(freqBins) length(timeBins) size(powMat,3)]);

% loop through freqBins
for f = 1:length(freqBins)
    %identify freq of interest
    [~,fIndEnd] = min(abs(freQ-freqBins(f)));
    if f == 1
        fIndBegin = 1;
    else
        [~,fIndBegin] = min(abs(freQ-freqBins(f-1)));
        fIndBegin = fIndBegin+1; %b/c the closest value to bin edge was already
                                 %used in the last iteration
    end
    fInd = fIndBegin:fIndEnd;
    %loop through timeBins
    for t = 1:length(timeBins)
        
        %identify time windows
        tInd = (timeBins(t)-timeWin+1):timeBins(t);
        
        %calculate mean power
        binPowMat(f,t,:) = nanmean(nanmean(powMat(fInd,tInd,:),1),2);
    end
    
end

%log 10 if flag is set
if logAfterMeanFlag ==1
    binPowMat = log10(binPowMat);
end
function [ev_out]=get_baseline_random(ev_in,secBetweenEv,secRandJitter,rSeed,recType)
%
% FUNCTION:
%   [ev_out]=getECoG_baseline_randTime(ev_in,secBetweenEv,secRandJitter,rSeed)
%
% DESCRIPTION:
%   gets randomly placed events in the eegfiles
%
% INPUTS:
%   ev_in............ events
%   secBetweenEv..... seconds between events (60)  
%   secRandJitter.... jitter around events times (10)
%   randSeed......... ranodm seed [optional]
%   recType..........'eeg' or other ('lfp' etc.) (optional; def = 'eeg')
%
% OUTPUTS:
%   ev_out........... baseline events
%
% NOTES:
%   (1) written by jfburke 10/2011 (john.fred.burke@gmail.com)
%   (2) AGR (ashwinramayya@gmail.com 04/2013).
%           -now works when session labels have more than one character
%           -additional input added: recording field (default 'eegfile'),
%           but can enter 'lfpfile' etc. if needed.

basEv = [];
counter=0;
%set the recType
if ~exist('recType','var')||isempty(recType)
    recType = 'eeg';
end

% set the random seed
if ~exist('rSeed','var')||isempty(rSeed)
  s = RandStream.create('mt19937ar','seed',sum(100*clock));
else
  s = RandStream.create('mt19937ar','seed',rSeed);
end
RandStream.setGlobalStream(s);

if isa(ev_in(1).session,'numeric')
    str2numFlag = 1;
    for i = 1:length(ev_in)
      ev_in(i).session = num2str(ev_in(i).session);
    end
else
    str2numFlag = 0;
end

uniSess    = unique({ev_in.session});
numSess    = length(uniSess);
thisSubj=unique({ev_in.subject});
useThisSubj=thisSubj{1};

for s=1:numSess  
  thisSess=uniSess(s);  
  if length({ev_in.session})~=length(ev_in);
    error('very bad');
  end  
  
  evThisSess = ev_in(strcmp({ev_in.session},thisSess));        
  eegFilesThisSess=unique({evThisSess.([recType 'file'])})';
  nFiles = length(eegFilesThisSess);
  
  % loop through events
  for f=1:nFiles
    thisEEGfile = eegFilesThisSess{f};
    if isempty(thisEEGfile);continue;end
  
    % get ino for this session
    [sampleFreq,nBytes,dataformat,gain] = GetRateAndFormat(thisEEGfile);
    sampleFreq = round(sampleFreq);
    if isempty(sampleFreq);error('NO PARAMS FILE');end
    
    % get the data type
    for k=1:1000
      filestem=sprintf('%s.%.3i',thisEEGfile, k);
      if ~exist(filestem,'file'); %try a different format
          filestem = sprintf('%s.%.i',thisEEGfile, k);
      end
      if ~exist(filestem,'file'); continue; end
      fileFeatures=dir(filestem);
      switch dataformat
       case {'int16','short'}
	numSamples=floor(fileFeatures.bytes/2);
       otherwise
	error('Are all dataformats int16 or short?')
      end
      break
    end
    
    intervalSam = floor(secBetweenEv*sampleFreq);    
    randomSam   = floor(secRandJitter*sampleFreq);    
    startSam    = intervalSam*2;
    endSam      = numSamples-intervalSam*2;  
    
    basSamplesThisFile_raw      = startSam:intervalSam:endSam;    
    nBasSamplesThisFile         = length(basSamplesThisFile_raw);  
    basSamplesThisFile_jittered = basSamplesThisFile_raw + ...	
	round((rand(1,nBasSamplesThisFile)-.5)*randomSam);    
    
    % see if you had any overflow    
    basSamplesThisFile_jittered(basSamplesThisFile_jittered<startSam)=startSam;    
    basSamplesThisFile_jittered(basSamplesThisFile_jittered>endSam)=endSam;    
    
    for k=1:nBasSamplesThisFile
      counter=counter+1;
      basEv(counter).subject=useThisSubj;
      basEv(counter).session=thisSess;
      basEv(counter).([recType 'file'])=thisEEGfile;
      basEv(counter).([recType 'offset'])=basSamplesThisFile_jittered(k);
    end    
    
  end
end
ev_out = basEv;

if str2numFlag == 1
   for i = 1:length(ev_out)
      ev_out(i).session = str2double(ev_out(i).session); 
   end
end    
% function [outStruct] = le_config_calcPow(configNum,task)
% %This function sets config for power analyses
% switch configNum
%     case 1
%         
%         
% end
% outStruct.configID = num2str(configNum);
% 
% 
% %constants:
% 
% %Freq stuff
% outStruct.freQ =  logspace(log10(2), log10(200),50);
% outStruct.freqBins = logspace(log10(2), log10(200),50);%[4 8 16 30 50 95];  
% outStruct.freqBinLbl={round((logspace(log10(2), log10(200),50)))};
% 
% %other param stuff
% outStruct.bipolFlag = 1;
% outStruct.zFlag = 1; 
% outStruct.filtFlag = 0; %60Hz filtFlag
% outStruct.logBefMeanFlag = 1; %
% outStruct.logAfterMeanFlag = 0;% these CANNOT both be set to 1 
% outStruct.bufferMS = 1000;
% outStruct.waveNum = 7;
% 
% %filter stuff
% outStruct.filtfreq = [];
% outStruct.filttype = [];
% outStruct.filtorder = [];
% 
% %'randBase' normalization params
% outStruct.baseSecBetEv = 10;
% outStruct.baseJitt = 5;
% outStruct.basePriorMS = 500;
% outStruct.basePostMS = 500;
% 
% % timing stuff compute offset and duration
% outStruct.timeWin = 500; 
% outStruct.timeStep = 100;
% if ~exist('task','var') || isempty(task)
%     outStruct.priorMS = 1000;
%     outStruct.postMS = 3000;
% else
%     switch task
%         case('pyFR')
%            outStruct.priorMS = 1000;
%            outStruct.postMS = 3000; 
%         case('prob_sel2')
%            outStruct.priorMS = 1000;
%            outStruct.postMS = 2000; 
%         case('motormap')
%            outStruct.priorMS = 1000;
%            outStruct.postMS = 5000;
%         case('CCDT')
%            outStruct.priorMS = 2000;
%            outStruct.postMS = 2500; 
%         case('trackball')
%            outStruct.priorMS = 1000;
%            outStruct.postMS = 3000; 
%         case{'tones3';'tones3.0'}
%             outStruct.zFlag = 0;
%             outStruct.timeWin = 100;
%             outStruct.timeStep = 20;
%             outStruct.priorMS = 0;
%             outStruct.postMS = 47600;
%         case{'LocGlobTones';'localGlobalTones'}
%             outStruct.timeWin = 400;
%             outStruct.timeStep = 200;
%             outStruct.priorMS = 0;
%             outStruct.postMS = 3700;
%             outStruct.freQ =  logspace(log10(0.3), log10(200),50);
%             outStruct.freqBins = logspace(log10(0.3), log10(200),50);
%             outStruct.baseSecBetEv = 8;
%             outStruct.baseJitt = 5;
%             outStruct.basePriorMS = 0;
%             outStruct.basePostMS = 3700;
%         case('oddball')
%             outStruct.priorMS = 500;
%             outStruct.postMS = 600; 
%     end
% end
% 
% outStruct.offsetMS = -outStruct.priorMS;
% outStruct.durationMS = outStruct.priorMS+outStruct.postMS;
% 
% % calculate time bins here for convenience
% tEnds = (outStruct.timeWin:outStruct.timeStep:outStruct.durationMS)+outStruct.offsetMS;
% tStarts = tEnds - outStruct.timeWin + 1;
% outStruct.timeBins = [tStarts' tEnds'];