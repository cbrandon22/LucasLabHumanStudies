function [dateStr,timeStr,ADchans,ADrefs]=parseThisSessionsLogFile(subj,rawDataPath,downsampleStruct)
%
% FUNCTION:
%  parseThisSessionsLogFile.m
%
% DESCRIPTION:
%
% INPUTS:
%  subj................ex: 'TJ030';
%  rawDataPath.........path from '/data/eeg/[subj]/raw' to data    
%  downsampleStruct....from nlx_downsample_config_file
%
% NOTES:
%  (1) written by jfburke 05/12 (john.fred.burke@gmail.com)
%  (2) This will only work with cheetah version 5.6
%
% CHEETAH DRS COMMANDS USAGE 
%  1. Channel Reference Select
%     This command will select one of the eight reference signals
%     to be used as the reference input for the specified output
%     channel (input board input signal).
%        DRS Command(bX crs ch# ref#)
%           bX  = board number
%           crs = Channel Reference Select option 
%           ch# = this is the channel number that will have its
%                 reference selected (values range from 0-31)
%           ref#=this is the reference signal that will be used for
%                the channel (values range from 0-7)
%  2. Reference Bus Source
%     This command is used to connect a headstage or global input
%     signal to one of the 8 reference bus signals.
%        DRS Command(bX rbs ref# chan#)
%           bX  = board number
%           rbs = Reference Bus Source option 
%           ref# = reference bus number valid values are 0-7
%           chan#=channel number that is to be used as the source
%                  0-31 select electrode inputs from the headstage
%                  32-35 select Ref1-Ref4 respectively
%                  36 selects animal ground
%                  37 selects Panel Ground
%                  38 selects the global ref input that corresponds
%                    to the ref# 
%


clear global
global CURRENTBUSINDS CURRENTBUSVALS ADCHANNELVECT ADCHANNELREFS ...
    NUMCHANNELPERBOARD DRS_BOARD_POSITION

logFileName        = downsampleStruct.logFileName;
NUMCHANNELPERBOARD = 32;
DRS_BOARD_POSITION = [];

% make the directories where the data live
dataDir = downsampleStruct.dataMotherDir;
subjDir = fullfile(dataDir,subj);
rawDir  = fullfile(subjDir,'raw',rawDataPath);
cfgFile = fullfile(rawDir,logFileName);
fid     = fopen(cfgFile,'r');
if fid==-1; 
  error(sprintf('\n\tCONFIG FILE ''%s'' doesn''t exist in %s',logFileName,rawDir));
end 

% these are the marker sthat let you know what line does what
startRecCommand     = 'AcquisitionControl::StartRecording()';
stopRecCommand      = 'AcquisitionControl::StopRecording()';
timeStampLineMarker = 'Time Opened';
drsCommandMarker    = 'HardwareSubSysDigitalLynxEdtInterface::SendDrsCommand()';
numPotChanMarker    = 'Licensed AD Channels';
addDrsBoardMarker   = 'Processing line: -AddDrsBoard';

% set up the bus lines, etc
numBusLines    = 8;
CURRENTBUSINDS = 0:numBusLines;
CURRENTBUSVALS = nan(1,numBusLines);

fprintf('START Parsing config file...\n')
while true
  N = fgetl(fid);
  if ~ischar(N);break;end

  %----------------------------------------------
  % get how many AD channels there potentially may be
  if ~isempty(regexp(N,numPotChanMarker))
    X = textscan(N,'%s%d','delimiter',':');
    
    % set up the bus lines, etc
    numPossibleChannels = double(X{2});
    ADCHANNELVECT       = 0:numPossibleChannels;
    ADCHANNELREFS       = nan(1,numPossibleChannels);
        
  %-----------------------------------------------
  % get order of DRS boards
  elseif ~isempty(regexp(N,addDrsBoardMarker));
    thisBoardNumber = str2double(N(end));
    if isnan(thisBoardNumber)
      error('bad drs board')
    end
    DRS_BOARD_POSITION = [DRS_BOARD_POSITION thisBoardNumber];
    
  %----------------------------------------------
  %        get time start/stop recording  
  elseif ~isempty(regexp(N,startRecCommand))||~isempty(regexp(N,stopRecCommand));;
    X = textscan(N,'%s%s%s%s%s%s','delimiter','-');
    if length(X{1})>1;error('bad config file read');end
    
    % get the time that you started the recording
    if ~isempty(regexp(N,startRecCommand))
      timeStr = getTimeStr_local(X{3}{1});
    end        
    fprintf('  %s at %s\n',X{6}{1},X{3}{1});        
  
    
  %----------------------------------------------
  %        get date of recording
  elseif ~isempty(regexp(N,timeStampLineMarker))
    X = textscan(N,'%s%s%s%s%s%s%s%s%s%s');
    dateStr = getDateStr_local(X{4}{1});
    fprintf('  Recording Date: %s\n',X{4}{1});
  
  %----------------------------------------------
  %  get what each channel is referenced to
  elseif ~isempty(regexp(N,drsCommandMarker))
    indOpenParantheses  = regexp(N,'(');
    indCloseParantheses = regexp(N,')');
    drsArguments        = N(indOpenParantheses(2)+1:end);
    X                   = textscan(drsArguments,'%s%s%d%d');
    thisBoard           = X{1}{1};
    thisDrsOption       = X{2}{1};
    switch upper(thisDrsOption)
     case {'RBS'}
      thisBusLine            = X{3};
      putThisAdChanOnTheLine = X{4};
      setBusLine(thisBoard,thisBusLine,putThisAdChanOnTheLine);
     case {'CRS'}
      thisAdChan            = X{3};
      thisAdChan_refBusLine = X{4};
      setCurrentRefs(thisBoard,thisAdChan,thisAdChan_refBusLine);
     case {'HSP'}
      %fprintf('  What does HSP mean?\n')
     otherwise
      error('bad drs option')
    end    
  else
    continue
  end
  
  
end
fclose(fid);
fprintf('END Parsing config file\n')

ADchans = ADCHANNELVECT;
ADrefs  = ADCHANNELREFS;

function setBusLine(board,bus,chan)
  global CURRENTBUSINDS CURRENTBUSVALS ADCHANNELVECT ADCHANNELREFS ...
      NUMCHANNELPERBOARD DRS_BOARD_POSITION 

  bNum = getBoardNumnerFromString(board);
  bInd = find(bNum);
  if isempty(bInd);
    error('cant find board');
  end
  ADstartChannel = (bInd-1)*NUMCHANNELPERBOARD;
  
  thisBusInd = find(bus==CURRENTBUSINDS);
  if isempty(thisBusInd)
    error('Bus index not found');
  end  
  CURRENTBUSVALS(thisBusInd)=ADstartChannel+chan;
  
function setCurrentRefs(board,chan,bus)
  global CURRENTBUSINDS CURRENTBUSVALS ADCHANNELVECT ADCHANNELREFS ...
      NUMCHANNELPERBOARD DRS_BOARD_POSITION 
  
  % get the board number and actual AD channe number (if this is
  % the second borad, then add that to the AD channel number)
  bNum = getBoardNumnerFromString(board);
  bInd = find(bNum);
  if isempty(bInd);
    error('cant find board');
  end
  ADstartChannel = (bInd-1)*NUMCHANNELPERBOARD;
  thisChannel    = chan + ADstartChannel;
    
  % get the value on the bus that this points to
  thisBusInd = find(bus==CURRENTBUSINDS);
  if isempty(thisBusInd)
    error('Bus index not found');
  end
  thisChannels_ref = CURRENTBUSVALS(thisBusInd);
    
  % store the value of this channels ref
  chanInd = find(ADCHANNELVECT==thisChannel);
  if isempty(chanInd)
    error('Channel index not found');
  end  
  ADCHANNELREFS(chanInd) = thisChannels_ref;
  
function bNum = getBoardNumnerFromString(bStr)
   bNum = str2double(bStr(end));
   if isnan(bNum)
     error('bad board number')
   end
  
function dateStr = getDateStr_local(dateStr_in);
  ind = regexp(dateStr_in,'/');
  m_str  = dateStr_in(1:ind(1)-1);
  d_str  = dateStr_in(ind(1)+1:ind(2)-1);
  y_str  = dateStr_in(ind(2)+3:end); 
  
  % make sure day has two digits
  d_str  = makeSureStrIsATwoDigitInteger_local(d_str);
  y_str  = makeSureStrIsATwoDigitInteger_local(y_str);  
  
  % get name of the month
  mm_int = str2double(m_str);
  if isnan(mm_int)
    error('The month is not an number.  Unexpected.')
  end    
  [~,mm]=month(sprintf('%02d',mm_int),'mm');  
  dateStr  = [d_str mm y_str];
  
function str_out = makeSureStrIsATwoDigitInteger_local(str_in)
  str_tmp = str2double(str_in);
  if int8(str_tmp)~=str_tmp; 
    error([str_in ' is not an integer less than 128'])
  end
  str_out = sprintf('%.2i',str_tmp);  
  
function timeStr = getTimeStr_local(timeStr_in);
  ind = regexp(timeStr_in,':');
  hour_str = timeStr_in(1:ind(1)-1);
  min_str  = timeStr_in(ind(1)+1:ind(2)-1);
  hour_str = makeSureStrIsATwoDigitInteger_local(hour_str);
  min_str  = makeSureStrIsATwoDigitInteger_local(min_str);     
  timeStr  = [hour_str min_str];

  