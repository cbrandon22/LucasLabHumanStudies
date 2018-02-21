function[subjPhase,config] = le_calcPhase(subj,sessList,elecNum,task,evLabel,configNum)
%% function[] = le_calcPhase(configNum,elecNum,subj,sess,task,evLabel)
%This function loads phase for a particular electrode or electrode pair
%associated with a particular subject's events from a particular task.
%Filters by session so that phase from relevant sessions are returned. 
%If phase is being computed for the first time, it saves it in the scratch
%directory defined in ec_dirs

%Inputs
%subj  = 'UP034'; 
%sessList = {'0'}; % (optional; if empty, it loops through all sessions seperately all events); 
                %Note: must be a cell array of strings
%elecNum = [1 2] (for bipolar contacts)
%task = 'pyFR'; %the label of the task dir in /data/events
%evLabel = 'math' (optional: default is 'events', but allows you to load
          %math or recog from pyFR if needed)

%Outputs
%subjPhase        concatonated phase for each events (freq x time x event); 
                %see config file to get time and freq windows
%Written by AGR 09-23-2013
%% set dirs and parse inputs
dirs = le_dirs;
if ~exist('evLabel','var') || isempty(evLabel)
    evLabel = 'events';
end
if ~exist('configNum','var') || isempty(configNum)
    configNum = 1;
end

%% load config
config = le_config_calcPhase(configNum,task);

%% go to pow directory
phaseDir = fullfile(dirs.scratch,'PHASE',num2str(config.configID),subj,[task '_' evLabel]);
cd_mkdir(phaseDir);
subjPhase = [];

%% load events
%load events
[events] = ec_loadEvents(subj,sessList,task,evLabel);
unique([events.session])
%% if ~exist sessList, create it here
if ~exist('sessList','var') || isempty(sessList)
    sessList = unique({events.session});
end

%% loop through sessList
for i = 1:length(sessList)
    sess = sessList{i};
    cd_mkdir(fullfile(phaseDir,sess))
    phaseFile = [num2str(elecNum(1)) '-' num2str(elecNum(2)) '.mat'];
    
    %try and load the phase file
    if exist(fullfile(phaseDir,sess,phaseFile),'file')
        load(phaseFile)
    else 
        %filter events
        sessEv = events(strcmp({events.session},sess));
        %calculate phase for this session and save
        phase = single(calcPhase(elecNum,sessEv,config));
        save(phaseFile,'phase')
    end
    
    %concatonate pow
    subjPhase = cat(3,subjPhase,phase);
    disp([num2str(elecNum) '-' sess]);
end

%% functions below
function[phase] = calcPhase(elecNum,events,config)
%this function calculates phase for a particular session
%gete_ms_wrapper
[EEG] = gete_ms_wrapper(elecNum,events, ...
               config.durationMS,config.offsetMS,config.bufferMS,...
               config.filtfreq, config.filttype, config.filtorder);
%eeg2pow 
[~,phase] = eeg2PowAndPhase(EEG,config.freQ,config.bufferMS,[],config.waveNum); 
phase = single(phase);

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
    for k=1:length(events)
         events(k).eegfile=regexprep(events(k).eegfile,'eeg.reref','eeg.noreref');
    end
    
    %Second, gete from both channels and subtract them from each other 
    [eeg1,resampfreq] = gete_ms(elecNum(1), events,durationMS,offsetMS,bufferMS,...
        filtfreq,filttype,filtorder,1000);
    [eeg2,resampfreq] = gete_ms(elecNum(2), events,durationMS,offsetMS,bufferMS,...
        filtfreq,filttype,filtorder,1000);
    
    eeg = eeg1 - eeg2;
end
function[powEEG,phaseEEG] = eeg2PowAndPhase(eeg,freQ,buffToRem,...
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

    
function [outStruct] = le_config_calcPhase(configNum,task)
%This function sets config for phase analyses
switch configNum
    case 1
        
        
end
outStruct.configID = num2str(configNum);


%constants:

%Freq stuff
outStruct.freQ =  logspace(log10(2), log10(200),50);
outStruct.freqBins = logspace(log10(2), log10(200),50);%[4 8 16 30 50 95];  
outStruct.freqBinLbl={round((logspace(log10(2), log10(200),50)))};

%other param stuff
outStruct.bipolFlag = 1;
outStruct.filtFlag = 0; %60Hz filtFlag
outStruct.bufferMS = 1000;
outStruct.waveNum = 7;

%filter stuff
outStruct.filtfreq = [];
outStruct.filttype = [];
outStruct.filtorder = [];

% timing stuff compute offset and duration
outStruct.timeWin = 500;
outStruct.timeStep = 100;
if ~exist('task','var') || isempty(task)
    outStruct.priorMS = 1000;
    outStruct.postMS = 3000;
else
    switch task
        case('pyFR')
           outStruct.priorMS = 1000;
           outStruct.postMS = 3000; 
        case('prob_sel2')
           outStruct.priorMS = 1000;
           outStruct.postMS = 2000; 
       case('motormap')
           outStruct.priorMS = 1000;
           outStruct.postMS = 5000; 
        case('trackball')
           outStruct.priorMS = 1000;
           outStruct.postMS = 3000; 
        case{'tones3';'tones3.0'}
            outStruct.zFlag = 0;
            outStruct.timeWin = 100;
            outStruct.timeStep = 20;
            outStruct.priorMS = 0;
            outStruct.postMS = 47600;
        case{'LocGlobTones';'localGlobalTones'}
            outStruct.timeWin = 400;
            outStruct.timeStep = 200;
            outStruct.priorMS = 0;
            outStruct.postMS = 3700;
            outStruct.freQ =  logspace(log10(0.3), log10(200),50);
            outStruct.freqBins = logspace(log10(0.3), log10(200),50);
            outStruct.baseSecBetEv = 8;
            outStruct.baseJitt = 5;
            outStruct.basePriorMS = 0;
            outStruct.basePostMS = 3700;
    end
end

outStruct.offsetMS = -outStruct.priorMS;
outStruct.durationMS = outStruct.priorMS+outStruct.postMS;


% calculate time bins here for convenience
tEnds = (outStruct.timeWin:outStruct.timeStep:outStruct.durationMS)+outStruct.offsetMS;
tStarts = tEnds - outStruct.timeWin + 1;
outStruct.timeBins = [tStarts' tEnds'];