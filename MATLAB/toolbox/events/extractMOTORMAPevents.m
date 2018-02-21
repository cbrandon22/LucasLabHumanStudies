function extractMOTORMAPevents(subject,expDir,session,forceSession)
%
% FUNCTION:
%   extractMOTORMAPevents(SUBJECT,EXPDIR,SESSION,FORCESESSION)
% 
% DESCRIPTION:
%   wrapper functions that makes and saves extractMOTORMAPevents events
%
% INPUTS:
%   SUBJECT.........'TJ051'
%   EXPDIR..........path to 'session_['SESSION']' directory.  Examples:
%                   '/data/eeg/TJ038_1/behavioral/motormap'
%                   '/data/eeg/TJ020/behavioral/motormap/TJ020c' 
%   SESSION.........0 = looks in 'session_0' in EXPDIR
%   FORCESESSION....1 = [optional] sets session to 1 (despite the
%                       fact that behavioral data are in session_0)
%                       Leave blank or empty if session number is
%                       same as SESSION
%
% OUTPUTS:
%   Makes and save three events:
%     (1) events.mat......... contains 'events'
%
% NOTES:
%   (1) written by jfburke on 2012-11-08 (john.fred.burke@gmail.com)
%

% check tp see if subject names match
fprintf('\n')
if isempty(regexp(expDir,subject))
  fprintf('  WARNNG: %s not found in %s\n',upper(subject),upper(expDir))
  fprintf('          you might be making an error.\n')
  fprintf('          please check this before making events. EXITING\n\n')
  fprintf('               !!! NO EVENTS SAVED !!! \n\n')
  return
end

% set defaults
if ~exist('forceSession','var')
  forceSession = [];
end

% get the directories
thisSessDirNAME = sprintf('session_%d',session);
thisSessDir     = fullfile(expDir,thisSessDirNAME);
evFile          = fullfile(thisSessDir,'events.mat');


% print opening line
fprintf('  Making MOTORMAP EVENTS for ')
fprintf('%s, session %d: \n',subject,session)

%--------------------------------------------
if ~exist(evFile,'file')
  events=extractMOTORMAPevents_local(subject,expDir,session,forceSession);
  if isempty(events)
    fprintf('BAD:  EMPTY EVENTS... exiting.\n\n')
    return
  end
  save(evFile,'events');
  clear events
  fprintf('DONE.\n')
else
  fprintf('SKIPPING (events exist).\n')
end

%----------------------------------------------
%fix the eeg log
fprintf('  Fixing EEG.EEGLOG file...')
eegLogFile_old = fullfile(thisSessDir,'eeg.eeglog');
eegLogFile_new = fullfile(thisSessDir,'eeg.eeglog.up');
if ~exist(eegLogFile_old,'file')
  fprintf('BAD:  NO EEG.EEGLOG FILE... exiting.\n\n')
  return
end
fixEEGLog(eegLogFile_old,eegLogFile_new);
pause(.2)
fprintf('done\n\n')

%<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>
%   ___orz__O7Z__]------------<>-**-<>------------[___orz__O7Z__
%<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>
function events=extractMOTORMAPevents_local(subject,expDir,session,forceSession);
  clear global
  global SUBJECT SESSION LIST MOVE events

  SUBJECT = subject;
  SESSION = session;
  LIST    = -999;
  MOVE    = 'X';
  
  thisSessDir = sprintf('session_%d',SESSION);
  sessFile    = fullfile(expDir,thisSessDir,'session.log');
  
  fid = fopen(sessFile,'r');
  if fid==-1
    fprintf('session %d..no session.log file found.\n',SESSION);   
    fprintf('EXITING\n\n');
    return
  end
  
  % you can change the session
  if exist('forceSESSION','var') && ~isempty(forceSESSION)
    SESSION=forceSESSION;
  end

  % get experimental variables
  numSess         = getExpInfo_local(expDir,'numSessions'); 
  numTrialPerCond = getExpInfo_local(expDir,'RepsPerCondition'); 
  motionCues      = getExpInfoStr_local(expDir,'ThingsToMove'); 
  numMotionCues   = length(motionCues);  
  numTrials       = numTrialPerCond*numMotionCues;
  
  evCounter         = 0;
  events            = [];
  sessLineCounter   = 0; 

  while true
    thisLine = fgetl(fid);
    if ~ischar(thisLine);return;end
    
    % get the third string before the underscore
    xTOT=textscan(thisLine,'%f%f%s','delimiter','\t');
    thisTYPE    = xTOT{3}{1};
    % based on the type write different fields for this event
    switch upper(thisTYPE)    
      %-------------------------------------------------------------------
     case {'B','E','SESS_END','TRIAL_END','WAIT','INSTRUCT','MOVE'}
      sessLineCounter = sessLineCounter + 1;
      x=textscan(thisLine,'%f%f%s','delimiter','\t');
      evCounter = evCounter + 1;
      mkNewEvent_local(evCounter,x{1},x{2});  
      appendNewEvent_local(evCounter,'type',x{3}{1});      
      %-------------------------------------------------------------------
     case 'SESS_START'
      sessLineCounter = sessLineCounter + 1;
      x=textscan(thisLine,'%f%f%s%d','delimiter','\t');
      thisSess  = x{4};
      if (thisSess-1)~=session
	error('sessions dont match');
      end
      if thisSess>numSess
	error('exceeds max session');
      end    
      evCounter = evCounter + 1;
      mkNewEvent_local(evCounter,x{1},x{2});  
      appendNewEvent_local(evCounter,'type',x{3}{1});       
      %-------------------------------------------------------------------
     case 'TRIAL_START'
      sessLineCounter = sessLineCounter + 1;
      x=textscan(thisLine,'%f%f%s%d%s','delimiter','\t');
      thisTRIAL = x{4};
      MOVE  = x{5}{1};
      if thisTRIAL>numTrials
        error('exceeds max trials')
      end
      LIST = thisTRIAL;
      evCounter = evCounter + 1;
      mkNewEvent_local(evCounter,x{1},x{2});  
      appendNewEvent_local(evCounter,'type',x{3}{1});
      %-------------------------------------------------------------------
     otherwise 
      error('type not recognized')
    end % end type switch
     
end % end fgetl

      
function mkNewEvent_local(evCounter,mstime,msoffset)
  global SUBJECT SESSION LIST MOVE events 

  events(evCounter).subject       = SUBJECT;
  events(evCounter).session       = SESSION;
  events(evCounter).trial         = LIST;
  events(evCounter).type          = 'X';
  events(evCounter).item          = MOVE;
  events(evCounter).mstime        = mstime;
  events(evCounter).msoffset      = msoffset;
  
function appendNewEvent_local(evCounter,varargin)
  global events
  nVar = length(varargin)/2;
  for v=1:nVar
    thisVarField      = varargin{2*(v-1)+1};
    thisVarData       = varargin{2*(v-1)+2};
    events(evCounter) = setfield(events(evCounter),thisVarField,thisVarData);
  end
  
function [out] = getExpInfo_local(expDir,str2get);
  fid_foo1 = fopen(fullfile(expDir,'config.py'),'r');
  while true
    thisLine = fgetl(fid_foo1);
    if ~ischar(thisLine);break;end
    if numel(thisLine)==0;continue;end
    if strcmp(thisLine(1),'#');continue;end
    possible_str=textscan(thisLine,'%s%f','Delimiter','=');
    X = regexprep(possible_str{1}{1},' ','');
    if strcmp(X,str2get)
      out=possible_str{2};
      break
    end
  end
  fclose (fid_foo1);
  
function [conds_all] = getExpInfoStr_local(expDir,str2get);
  fid_foo1 = fopen(fullfile(expDir,'config.py'),'r');
  while true
    thisLine = fgetl(fid_foo1);
    if ~ischar(thisLine);break;end
    if numel(thisLine)==0;continue;end
    if strcmp(thisLine(1),'#');continue;end
    possible_str=textscan(thisLine,'%s%s','Delimiter','=');
    X = regexprep(possible_str{1}{1},' ','');
    if strcmp(X,str2get)
      out=possible_str{2};
      break
    end
  end
  fclose (fid_foo1);
  
  switch upper(str2get)    
   case 'THINGSTOMOVE'
    conds_all = [];
    thisStr = out{1};
    counter=0;
    firstApostrophe=false;
    for k=1:length(thisStr)
      if strcmp(thisStr(k),'"')&~firstApostrophe
	counter=counter+1;
	firstApostrophe=true;
	conds_tmp = cell(1,1);
	continue
      elseif strcmp(thisStr(k),'"')&firstApostrophe
	firstApostrophe=false;
	conds_all = cat(1,conds_all,conds_tmp);
	continue
      else
	if firstApostrophe
	  conds_tmp{1}(end+1)=thisStr(k);
	end
      end
    end
    
   otherwise
    error('bad variable')
  end
  