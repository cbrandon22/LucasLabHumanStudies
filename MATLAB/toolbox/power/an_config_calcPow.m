function [outStruct] = an_config_calcPow(configNum)
%This function sets config for power analyses
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
outStruct.bipolFlag = 0;
outStruct.zFlag = 1; 
outStruct.filtFlag = 0; %60Hz filtFlag
outStruct.logBefMeanFlag = 1; %
outStruct.logAfterMeanFlag = 0;% these CANNOT both be set to 1 
outStruct.bufferMS = 1000;
outStruct.waveNum = 7;
outStruct.nlxFlag = 0;

%filter stuff
outStruct.filtfreq = [];
outStruct.filttype = [];
outStruct.filtorder = [];

%'randBase' normalization params
outStruct.baseSecBetEv = 10;
outStruct.baseJitt = 5;
outStruct.basePriorMS = 500;
outStruct.basePostMS = 500;

% timing stuff compute offset and duration
outStruct.timeWin = 500; 
outStruct.timeStep = 100;
if ~exist('task','var') || isempty(task)
    outStruct.priorMS = 1000;
    outStruct.postMS = 3000;
else
    switch task
        case('oddball')
           outStruct.priorMS = 1000;
           outStruct.postMS = 5000;
           outStruct.nlxFlag = 1;
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