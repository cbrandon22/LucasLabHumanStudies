function extractGONOGOevents(subject,expDir,session,forceSession)
%
% FUNCTION:
%   extractMOTORMAPevents(SUBJECT,EXPDIR,SESSION,FORCESESSION)
% 
% DESCRIPTION:
%   wrapper functions that makes and saves extractGONOGOevents events
%
% INPUTS:
%   SUBJECT.........'TJ051'
%   EXPDIR..........path to 'session_['SESSION']' directory.  Examples:
%                   'gonogo/data/SUBJ/'
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
%   (2) GONOGO version by cbrandon on 10/13/2017 (cbrandon22@gmail.com)

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
fprintf('  Making GONOGO EVENTS for ')
fprintf('%s, session %d: \n',subject,session)

%--------------------------------------------
if ~exist(evFile,'file')
  events=extractODDBALLevents_local(subject,expDir,session,forceSession);
  if isempty(events)
    fprintf('File not found:  Dir should contain session folder.\n\n')
    return
  end
  save(evFile,'events');
  clear events
  fprintf('DONE.\n')
else
  fprintf('SKIPPING (events exist).\n')
end

%<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>
%   ___orz__O7Z__]------------<>-**-<>------------[___orz__O7Z__
%<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>
function events=extractODDBALLevents_local(subject,expDir,session,forceSession);
  clear global
  global SUBJECT SESSION LIST RESP events

  SUBJECT = subject;
  SESSION = session;
  LIST    = -999;
  RESP    = -999;
  
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
  numBlocks       = getExpInfo_local(expDir,'blocksPerSession'); 
  trialsPerBlock  = getExpInfo_local(expDir,'trialsPerBlock');
  flipTones       = getExpInfo_local(expDir,'flipTones');
  numTrials       = numBlocks*trialsPerBlock;
  
  evCounter         = 0;
  events            = [];
  sessLineCounter   = 0;
  passiveCount      = 0;

  while true
    thisLine = fgetl(fid);
    if ~ischar(thisLine);return;end
    
    % get the third string before the underscore
    xTOT=textscan(thisLine,'%f%f%s','delimiter','\t');
    thisTYPE    = xTOT{3}{1};
    % based on the type write different fields for this event
    if strcmp(thisTYPE,'PASSIVE_START')
        passiveCount = passiveCount+1;
        rblock = 0;
        targTone = 1;
        if strcmp(flipTones,'True') && rem(passiveCount,2)==0
            targTone = 2;
        end
    elseif strcmp(thisTYPE,'RESPONSE_START')
        rblock = 1;
        targTone = 0;
    end
    switch upper(thisTYPE)    
      %-------------------------------------------------------------------
     case {'B','E','SESS_END','START_INSTRUCTIONS','DONE_INSTRUCTIONS','INSTRUCTION_AUDIO','PASSIVE_START','END_SOUND','RESPONSE_START','PASSIVE_END','RESPONSE_END'}
      sessLineCounter = sessLineCounter + 1;
      x=textscan(thisLine,'%f%f%s%d','delimiter','\t');
      RESP = -999;
      evCounter = evCounter + 1;
      mkNewEvent_local(evCounter,x{1},x{2});  
      appendNewEvent_local(evCounter,'type',x{3}{1});      
      %-------------------------------------------------------------------
     case 'SESS_START'
      sessLineCounter = sessLineCounter + 1;
      x=textscan(thisLine,'%f%f%s%d','delimiter','\t');
      thisSess  = x{4};
      RESP = -999;
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
     case {'BACKGROUND','TARGET'}
      sessLineCounter = sessLineCounter + 1;
      x=textscan(thisLine,'%f%f%s%d','delimiter','\t');
      thisTRIAL = x{4};
      RESP = -999;
      if thisTRIAL>numTrials
        error('exceeds max trials')
      end
      LIST = thisTRIAL;
      evCounter = evCounter + 1;
      mkNewEvent_local(evCounter,x{1},x{2});
      if rblock && strcmp(thisTYPE,'BACKGROUND')
          x{3}{1} = 'NO_GO';
      elseif rblock && strcmp(thisTYPE,'TARGET')
          x{3}{1} = 'GO';
      elseif targTone ==1 && strcmp(thisTYPE,'BACKGROUND')
          x{3}{1} = 'BACKGROUNDLF';
      elseif targTone ==1 && strcmp(thisTYPE,'TARGET')
          x{3}{1} = 'TARGETHF';
      elseif targTone ==2 && strcmp(thisTYPE,'BACKGROUND')
          x{3}{1} = 'BACKGROUNDHF';
      elseif targTone ==2 && strcmp(thisTYPE,'TARGET')
          x{3}{1} = 'TARGETLF';
      end
      appendNewEvent_local(evCounter,'type',x{3}{1});
      %-------------------------------------------------------------------
     case 'RESPONSE'
      sessLineCounter = sessLineCounter + 1;
      x=textscan(thisLine,'%f%f%s%d%s%s','delimiter','\t');
      %check response T\F
      respBool = regexp(x{6},'=','split');
      respBool = respBool{1}{2};
      if strcmp(respBool,'False')
          RESP = 0;
      elseif strcmp(respBool,'True')
          RESP = 1;
      else
          error(['unrecognized response: ' respBool]);
      end
      evCounter = evCounter + 1;
      mkNewEvent_local(evCounter,x{1},x{2});
      appendNewEvent_local(evCounter,'type',x{3}{1});
      %-------------------------------------------------------------------
     otherwise 
      error('type not recognized')
    end % end type switch
     
end % end fgetl

      
function mkNewEvent_local(evCounter,mstime,msoffset)
  global SUBJECT SESSION LIST RESP events 

  events(evCounter).subject       = SUBJECT;
  events(evCounter).session       = SESSION;
  events(evCounter).trial         = LIST;
  events(evCounter).type          = 'X';
  events(evCounter).response      = RESP;
  events(evCounter).mstime        = mstime;
  events(evCounter).msoffset      = msoffset;
  if RESP ~= -999
      events(evCounter-1).response = RESP;
      events(evCounter-2).response = RESP;
  end
  
function appendNewEvent_local(evCounter,varargin)
  global events
  nVar = length(varargin)/2;
  for v=1:nVar
    thisVarField      = varargin{2*(v-1)+1};
    thisVarData       = varargin{2*(v-1)+2};
    events(evCounter) = setfield(events(evCounter),thisVarField,thisVarData);
  end
  
function [out] = getExpInfo_local(expDir,str2get)
  fid_foo1 = fopen(fullfile(expDir,'config.py'),'r');
  while true
      strFlag=0;
    thisLine = fgetl(fid_foo1);
    if ~ischar(thisLine);break;end
    if numel(thisLine)==0;continue;end
    if strcmp(thisLine(1),'#');continue;end
    possible_str=textscan(thisLine,'%s%f','Delimiter','=');
    if isempty(possible_str{2}) % return as string
        possible_str=textscan(thisLine,'%s%s','Delimiter','=');
        strFlag=1;
    end
    X = regexprep(possible_str{1}{1},' ','');
    if strcmp(X,str2get)
        if strFlag
            out=possible_str{2}{1};
        else
            out=possible_str{2};
        end
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