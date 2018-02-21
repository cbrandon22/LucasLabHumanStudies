function [fileExt] = an_downsampleAndSaveNlxData(subj,rawDataPath,logpeg)
% this function processes the raw nlx files and saves out a downsampled 
% lfp file for each recorded channel in the lfp.noreref folder in each
% subject's directory
% FUNCTION:
%  an_downsampleAndSaveNlxData.m
%
% DESCRIPTION:
%
% INPUTS:
%  subj= ex: 'TJ030';
%  rawDataPath= path from data directory '/data/eeg/[subj]/raw' to data    
%
% OUTPUTS:
% fileExt....filename of lfp.noreref file
% NOTES:
%  (1) written by jfburke 05/12 (john.fred.burke@gmail.com)
% 
% To-Do:
%   (1) get the reference up there
%   (2) get the following fileds from stat_ncs.m 
%       i.  InputInverted
%       ii. DSPLowCutFilterEnabled
%       ii. DSPHighCutFilterEnabled
%

% initialize
if logpeg == 1
    an_nlx_downsample_config_file2
else
    an_nlx_downsample_config_file % yields downsampleStruct
end

% make the directories where the data live
dataDir = downsampleStruct.dataMotherDir;
subjDir = fullfile(dataDir,subj);
rawDir  = fullfile(subjDir,'raw',rawDataPath);
lfpDir  = fullfile(subjDir,'lfp.noreref');
if ~exist(lfpDir);mkdir(lfpDir);end

% get all ncs files with data
[elecNames, evFileName] = getAllCscFilesWithActualData(rawDir,downsampleStruct);
nElecs                 = length(elecNames);
fprintf('\n')
fprintf('%s: \n',subj)
fprintf('%s\n',rawDir);pause(.2);% because I like it, trick
fprintf('Found %d files\n',nElecs);pause(.2);
fprintf('Extracting all channels....\n\n\n')

% parse log file for info
[dateStr,timeStr,AD_Channels_all,AD_references_all]=...
    parseThisSessionsLogFile(subj,rawDataPath,downsampleStruct);
fprintf('\n')

% define globals for use in this function
clear global
global DOWNSAMPLE_TSDIFF DOWNSAMPLE_RATE LOWPASSFREQ FIRSTSAMPLETIME ...
    LASTSAMPLETIME ... 
    NUMDOWNSAMPLES DATA_TYPE ADCHANNEL INPUTRANGE ...
    DSPLOWCUTFREQUENCY DSPHIGHCUTFREQUENCY DRS_REFERENCE_AD ...
    ORIG_SAMPLERATE

% make the file extension
fileExt = sprintf('%s_%s_%s',subj,dateStr,timeStr);
fprintf('File Extension: %s\n',fileExt);pause(.2);

% open and initialize the jacksheet
jackSheetFileName = sprintf('%s.JacksheetAndParams.txt',fileExt);
jackSheetFile     = fullfile(lfpDir,jackSheetFileName);
if exist(jackSheetFile,'file');
  fprintf('\n\n\t!!!EXITING!!!: %s exists already in %s\n\n',jackSheetFileName,lfpDir)
  return
end
fid=fopen(jackSheetFile,'w','l');
if fid==-1;
  error(sprintf('Cannot open %s in %s',jackSheetFileName,lfpDir));
end
initializeJackksheet_local(fid);

% loop through each channel
fprintf('\n\n')
for k=1:length(elecNames)

  % see if it is already done
  cscNum = str2double(regexp(elecNames{k},'\d*','match'));
  chanfileName = sprintf('%s.%03i', fileExt,cscNum);
  chanfile     = fullfile(lfpDir,chanfileName);
  if exist(chanfile,'file')
    fprintf('  %s already exists... SKIPPING\n',chanfileName)
    continue
  end
  
  % load the data
  thisChannelFile = fullfile(rawDir,elecNames{k});
  if ~exist(thisChannelFile,'file');
    error('This makes no sense and should never happen')
  end  
  fprintf('CHANNEL %d: %s\n',k,elecNames{k})
  fprintf('---------------------\n')
  [data,info,TSblock,TSsamples]=load_ncs2(thisChannelFile);
  ADCHANNEL  = info.ADChannel;
  INPUTRANGE = info.InputRange;
  adChanInd = find(AD_Channels_all==ADCHANNEL);
  if isempty(adChanInd); error(sprintf('cannot find AD %d',ADCHANNEL));end
  DRS_REFERENCE_AD = AD_references_all(adChanInd);
  
  % get the number of samples to skip (to get closest to the desird
  % sampling rate without going under)
  if strcmp(downsampleStruct.targetSampleRate,'native')
      downsampleStruct.targetSampleRate = floor(info.actualSampleRate);
  end
  timeInSecBetweenActualSamples = (info.sampleTSDiff./1e6);  
  decimateFactor = floor(1./(downsampleStruct.targetSampleRate*timeInSecBetweenActualSamples));
  timeInSecBetweenDownsamples = decimateFactor*timeInSecBetweenActualSamples;
  DOWNSAMPLE_TSDIFF           = timeInSecBetweenDownsamples;
  DOWNSAMPLE_RATE             = 1./timeInSecBetweenDownsamples;
  DSPLOWCUTFREQUENCY          = info.DspLowCutFrequency;         
  DSPHIGHCUTFREQUENCY         = info.DspHighCutFrequency;
  ORIG_SAMPLERATE             = info.actualSampleRate;
    
  % low-pass filter the data to prevent aliasing
  LOWPASSFREQ = floor(downsampleStruct.engineersNyquist*DOWNSAMPLE_RATE);
  
  % convert data to doubles before filtering....
  fprintf('Converting data to doubles for filtering....')
  data = double(data);
  fprintf('done\n')
  
  % filter the data
  fprintf('Filtering above %d Hz....',LOWPASSFREQ)
  switch upper(downsampleStruct.lowPassType)
   case 'BUTTER'   
    
    T_min   = 20; % minutes of ECoG to run th filter
    T_sec   = T_min*60;
    data_LP = buttfilt(data,LOWPASSFREQ,info.actualSampleRate,'low',downsampleStruct.lowPassOrder);
    
    
   otherwise
    error('filter type not supported')
  end
  fprintf('done\n')

  % Perform the decimation
  fprintf('Decimating every %d samples.\n',decimateFactor)  
  indToKeep                        = false(1,length(data_LP));
  indToKeep(1:decimateFactor:end)  = true;
  downSampledData                  = data_LP(indToKeep);
  downSampledTS                    = TSsamples(indToKeep); 
  FIRSTSAMPLETIME = downSampledTS(1);
  LASTSAMPLETIME  = downSampledTS(end);
  NUMDOWNSAMPLES  = length(downSampledData);
  fprintf('New downsampling rate is %0.1f Hz\n',DOWNSAMPLE_RATE);
    
  % convert the downsampled data to int16's (the original data format)  
  fprintf('Converting filtered data back to int16....')
  DATA_TYPE             = 'int16';
  downSampledData_int16 = int16(downSampledData);
  fprintf('done\n')  
    
  % write the jacksheet line
  writeJackksheetLine_local(fid,cscNum,elecNames{k});
  
  % finally write the file
  fprintf('writing to %s.....',chanfile);
  outfid=fopen(chanfile,'w','l');
  c=fwrite(outfid,downSampledData_int16,'int16');
  fclose(outfid);
  fprintf('done\n')  
  
  clear data data_LP
  fprintf('\n\n\n')
    
end

% first open the sync pulse file
pulseFileName = sprintf('%s.sync.txt',fileExt);
pulseFile     = fullfile(lfpDir,pulseFileName);
if exist(pulseFile,'file');
  fprintf('\n\n\t!!!EXITING!!!: %s exists already in %s\n\n',pulseFileName,lfpDir)
  return
end

% now get and write the pulses
evFile = fullfile(rawDir,evFileName);
if ~exist(evFile,'file');
  fprintf('\n\n\tTHIS SHOULD NEVER HAPPEN: no %s in %s\n\n',evFileName,rawDir);
  return
end

[upstrokesTS]       = getNlxEventTimes(evFile);
%firstSampleTime_ALL = getNlxJackSheetInfo(lfpDir,fileExt,'FIRSTSAMPLETIME','all');
%upstrokesTS=upstrokesTS-firstSampleTime_ALL;
%eventsSamples=round(upstrokesTS/DOWNSAMPLE_TSDIFF)+1;
eventsSamples=upstrokesTS;
pulseFid=fopen(pulseFile,'w','l');
fprintf(pulseFid,'%f\n',eventsSamples);
fclose(pulseFid);

% write params.txt file (to make it compatible with gete)
cd(lfpDir)
paramsfid = fopen('params.txt','w','l');
fprintf(paramsfid,'samplerate\t%d\n',DOWNSAMPLE_RATE);
fprintf(paramsfid,'dataformat\t''%s''\n',DATA_TYPE);
fclose(paramsfid)
clear global

%-------------------------------------------------------------------------
function [allElecNames evFile] = getAllCscFilesWithActualData(rawSubjDir,downsampleStruct)
d        = dir(rawSubjDir);
thisExt  = downsampleStruct.fileExt;

% get events file
evFile = getEventsFile_local(d,downsampleStruct); 

% get all of the raw data files
rawDataStruct = d(~cellfun(@isempty,regexp({d.name},[thisExt '$'])));
rawDataFiles  = {rawDataStruct.name};
rawDataSizes  = [rawDataStruct.bytes];
if sum([rawDataStruct.isdir])>0;error('data ~= dir');end
gdInd         = rawDataSizes>downsampleStruct.minBytesInFile;
allElecNames  = rawDataFiles(gdInd);
allElecSizes  = rawDataSizes(gdInd);

errorCheck = std(allElecSizes)./mean(allElecSizes);
if errorCheck>0
  fprintf('  WARNING!! FILES ARE NOT ALL THE SAME SIZE\n')
  fprintf('            percent variation = %2.2f%%\n',errorCheck)
end

%-------------------------------------------------------------------------
function evFile_out=getEventsFile_local(d,exStruct)
  evFileInfo = d(strcmp({d.name},exStruct.evFile));
  if isempty(evFileInfo)
    fprintf(' WARNING!! %s FILE NOT FOUND\n')
    evFileInfo=[];return
  end
  evFile_out  = evFileInfo.name;
  evFile_size = evFileInfo.bytes;
  if evFile_size<=exStruct.minBytesInFile
    fprintf(' WARNING!! EVENTS LIKELY NOT RECORDED\n')
  end
    
%-------------------------------------------------------------------------  
function initializeJackksheet_local(fid);
global DOWNSAMPLE_TSDIFF DOWNSAMPLE_RATE LOWPASSFREQ FIRSTSAMPLETIME LASTSAMPLETIME ...
    NUMDOWNSAMPLES DATA_TYPE ADCHANNEL INPUTRANGE ...
    DSPLOWCUTFREQUENCY DSPHIGHCUTFREQUENCY DRS_REFERENCE_AD ...
    ORIG_SAMPLERATE
fprintf(fid,'channel_num\t');
fprintf(fid,'channel_name\t');
fprintf(fid,'ADCHANNEL\t');
fprintf(fid,'DRS_REFERENCE_AD\t');
fprintf(fid,'INPUTRANGE\t');
fprintf(fid,'DSPLOWCUTFREQUENCY\t');
fprintf(fid,'DSPHIGHCUTFREQUENCY\t');
fprintf(fid,'ORIG_SAMPLERATE\t');
fprintf(fid,'DOWNSAMPLE_TSDIFF\t');
fprintf(fid,'DOWNSAMPLE_RATE\t');
fprintf(fid,'LOWPASSFREQ\t');
fprintf(fid,'FIRSTSAMPLETIME\t');
fprintf(fid,'LASTSAMPLETIME\t');
fprintf(fid,'NUMDOWNSAMPLES\t');
fprintf(fid,'DATA_TYPE\t');
fprintf(fid,'\n');

%-------------------------------------------------------------------------  
function writeJackksheetLine_local(fid,channel_num,channel_name);
global DOWNSAMPLE_TSDIFF DOWNSAMPLE_RATE LOWPASSFREQ FIRSTSAMPLETIME LASTSAMPLETIME ...
    NUMDOWNSAMPLES DATA_TYPE ADCHANNEL INPUTRANGE ...
    DSPLOWCUTFREQUENCY DSPHIGHCUTFREQUENCY DRS_REFERENCE_AD ...
    ORIG_SAMPLERATE
fprintf(fid,'%d\t',channel_num);
fprintf(fid,'%s\t',channel_name);
fprintf(fid,'%d\t',ADCHANNEL);
fprintf(fid,'%d\t',DRS_REFERENCE_AD);
fprintf(fid,'%d\t',INPUTRANGE);
fprintf(fid,'%f\t',DSPLOWCUTFREQUENCY);
fprintf(fid,'%f\t',DSPHIGHCUTFREQUENCY);
fprintf(fid,'%f\t',ORIG_SAMPLERATE);
fprintf(fid,'%f\t',DOWNSAMPLE_TSDIFF);
fprintf(fid,'%f\t',DOWNSAMPLE_RATE);
fprintf(fid,'%d\t',LOWPASSFREQ);
fprintf(fid,'%f\t',FIRSTSAMPLETIME);
fprintf(fid,'%f\t',LASTSAMPLETIME);
fprintf(fid,'%d\t',NUMDOWNSAMPLES);
fprintf(fid,'%s\t',DATA_TYPE);
fprintf(fid,'\n');