function an_alignNeuralynxFile(fullFileExt,sessionDir,evCellArray,timesfile)
%
% FUNCTION:
%  alignNeuralynxFile.m
%
% DESCRIPTION:
% This function aligns 
%
% INPUTS:
%  fullFileExt = this is the fullfile name for the lfp file;
%  sessionDir = location of behavioral pulse files (eeg.eeglog.up)
%  evCellArray = cell array of behavioral events files (events.mat)
%  unit 
%
% NOTES:
%  (1) written by jfburke 05/12 (john.fred.burke@gmail.com)
%  (2) updated by AGR 04/13 
%              -edited preamble
%              -events now also add lfpoffset(ms), timesfile, and
%              timesoffset(ms) fields
% To-Do:
% ( ) Write Params file w/ sample rate

clear global
global FULLFILEEXT TIMESFILE
FULLFILEEXT = fullFileExt; 

if ~exist('timesfile','var') || isempty(timesfile)
    TIMESFILE = '';
else
    TIMESFILE = timesfile;
end
% 1. get the received eeg pulses
eegPulseFile = [fullFileExt '.sync.txt'];
if ~exist(eegPulseFile,'file')
  error(sprintf('%s does not exist\n',eegPulseFile))
end
[~,pulse_str] = system(['cut -f 1 ' eegPulseFile]);
eegPulse_ms = strread(pulse_str,'%f')./1e3; %sample numbers of received pulses

% 2. get the sent beahvioral pulses
behPulseFileName = 'eeg.eeglog.up';
behPulseFile     = fullfile(sessionDir,behPulseFileName);
behPulse_ms      = textread(behPulseFile,'%n%*[^\n]','delimiter','\t');

% 3. get the samplerate
[lfpDir,fileExt] = fileparts(fullFileExt);
samplerate = getNlxJackSheetInfo(lfpDir,fileExt,'ORIG_SAMPLERATE','all');
lfp_samplerate = getNlxJackSheetInfo(lfpDir,fileExt,'DOWNSAMPLE_RATE','all');

first_sample_time = getNlxJackSheetInfo(lfpDir,fileExt,'FIRSTSAMPLETIME','all');

% 4. align the pulses
threshMS  = 10;
window    = 500;
doplot    = true;
pulseIsMS = true; 
fprintf('\n')
fprintf('SYNCING\n')
fprintf('-------\n')
fprintf('%s --AND-- %s:\n',fileExt,sessionDir);
[behAligned_ms,nlxAligned_ms] = pulsealign(behPulse_ms,eegPulse_ms,samplerate,...
				 threshMS,window,doplot,pulseIsMS);
fprintf('Alignment: %d matches of %d recorded pulses',...
	     length(behAligned_ms),length(eegPulse_ms)); 
fprintf('\n\n')

% run logalign with the beh_ms and eeg_offset info
[a,err1]=logalign_local(behAligned_ms,nlxAligned_ms,evCellArray,samplerate,lfp_samplerate,first_sample_time);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [meddev,rawDevs]=logalign_local(beh_ms,nlx_ms,evCellArray,samplerate,lfp_samplerate,first_sample_time)
  %LOGALIGN_LOCAL - Fit line to EEG sync and add files and offsets to logs.
  %
  % You provide matched EEG(offsets) and behavioral(ms) sync pulses
  % and this function regresses a line for conversion from ms to eeg
  % offsets.  Then it goes through each log file and appends the
  % correct eeg file and offset into that file for each one.  You
  % must specify a sample eeg file for each set of eeg offsets so
  % that it can calculate the duration of each file.
  %

  global FULLFILEEXT TIMESFILE
  doplot = false;
  
  % allocate for stuff
  b = zeros(2,length(nlx_ms));
  eeg_start_ms = zeros(length(nlx_ms),1);
  eeg_stop_ms = zeros(length(nlx_ms),1);
  
  if doplot
    plot(beh_ms,nlx_ms,'.k','MarkerSize',15)
    xlabel('behavioral pulse times (msec)')
    ylabel('neuralynx pulse times (msec)')
  end
  
  % get slope and offset for each eeg file
  bfix = beh_ms(1);
  [b,bint,r,rin,stats] = regress(nlx_ms,[ones(length(beh_ms),1) beh_ms-bfix]);        
  b(1) = b(1) - bfix*b(2);
  
  if doplot
    hold on
    dummyX = linspace(beh_ms(1),beh_ms(end),1000)';
    dummyY = [ones(length(dummyX),1) dummyX]*b;
    plot(dummyX,dummyY,'--b','LineWidth',2)
    hold off
  end  
  
  % calc max deviation
  act=[ones(length(beh_ms),1) beh_ms]*b;
  maxdev = max(abs(act - nlx_ms));
  meddev = median(abs(act - nlx_ms));  
  rawDevs=act - nlx_ms;
      
  % report stats
  fprintf('\tMax. Dev. = %f ms\n', maxdev);
  fprintf('\tMedian. Dev. = %f ms\n', meddev);  
  fprintf('\t95th pctile. = %f ms\n', prctile(rawDevs,95)); 
  fprintf('\t99th pctile. = %f ms\n', prctile(rawDevs,99)); 
  fprintf('\tR^2 = %f\n', stats(1));
  fprintf('\tSlope = %f\n', b(2));
  fprintf('\tPulse range = %.3f minutes\n',range(beh_ms)/1000/60);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % now get all the event times (in beh_ms) and find where they lie
  % on the regression line (which will give nlx_ms)   
  numEventsToAlign = length(evCellArray);
  %cols             = myJet(numEventsToAlign);
  evTypes           = {}; 
  
  % loop through all events  
  for e=1:numEventsToAlign
    thisEv = evCellArray{e};
    if ~exist(thisEv,'file')
      fprintf('%s doesn''t exist\n',thisEv);
      continue
    end
    [~,fooTMP]=fileparts(thisEv);
    evTypes{e}=regexprep(fooTMP,'_',' ');
    
    % load the events
    events    = load(thisEv); 
    events    = events.events;
    behEvents_times = [events.mstime]'; 
    
    % find the nlx times
    nlxEvents_times = [ones(length(behEvents_times),1) behEvents_times]*b;    
    if doplot
      hold on
      h(e)=plot(behEvents_times,nlxEvents_times,...
		'x','Color',cols(e,:),'MarkerSize',10);
      hold off
    end
 
    % now fill in the nlx times for each event    
    numEvents = length(events);
    if ~isequal(numEvents,length(nlxEvents_times))
      error('should never happen')
    end
    for k=1:numEvents
      thisNlxEvTime = nlxEvents_times(k);      
      events(k).lfpfile = FULLFILEEXT;    
      events(k).lfpmstime = thisNlxEvTime;
      events(k).lfpoffset = round((thisNlxEvTime-(first_sample_time/1000))*(lfp_samplerate/1000)); %in lfp samples 
      events(k).timesfile = TIMESFILE;
      events(k).timesoffset = thisNlxEvTime;
    end
    if ~isequal([events.lfpmstime],nlxEvents_times')
      error('should never happen')
    end
    
    % now save the events and keep a copy of the old ones
    % 1. first create a copy
    evFile   = thisEv;
    cpEvFile = [thisEv '.lfp.old'];
    if exist(cpEvFile,'file')
      fprintf('%s is already aligned.... SKIPPING\n',evFile)
      continue
    end
    unixStr  = sprintf('mv %s %s',evFile,cpEvFile);
    unixBool = system(unixStr);
    if unixBool~=0
      error('Bad ''system'' call: Unable to make a copy of the events')
    end
    
    % 2. now save the newly aligned events
    save(evFile,'events');
    clear events evFile cpEvFile unixStr unixBool    
  end

  if doplot
    % make a legend
    legend(h,evTypes,'Location','NorthWest');
  end
  
function cols = myJet(n)
%get n random colors from the current colormap
cmap = colormap;
cols = cmap(randi([2],[1 n])-1,:);  
