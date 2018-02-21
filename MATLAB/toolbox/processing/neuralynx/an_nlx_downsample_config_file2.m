%This script is the config file that is used by an_downsampleAndSaveNlxData

[~,dataDir] = an_subjInfo2;
downsampleStruct.dataMotherDir     = dataDir; clear dataDir
%downsampleStruct.dataMotherDir     = '/data3/scratch/ramayya/anesthesia/';
downsampleStruct.logFileName       = 'PegasusLogFile.txt'; 
downsampleStruct.targetSampleRate  = 'native'; % 2000; 32768;
downsampleStruct.recType           = 'pegasus'; % cheetah
downsampleStruct.evFile            = 'Events.nev';
downsampleStruct.fileExt           = 'ncs';
downsampleStruct.minBytesInFile    = 20000;
downsampleStruct.engineersNyquist  = 2/5;
downsampleStruct.lowPassType       = 'BUTTER';
downsampleStruct.lowPassOrder      = 4;