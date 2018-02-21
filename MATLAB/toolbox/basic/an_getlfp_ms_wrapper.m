function[lfp] = an_getlfp_ms_wrapper(elecNum,events, ...
               durationMS,offsetMS,bufferMS,filtfreq,filttype,filtorder)
%This wrapper returns raw lfp traces for a set of events. 
%Inputs
%elecNum        %electrode number to load. Can be a single channel (eg,[1]) or a set
                %of channels (for bipolar, eg,[1 2]). see getBipolarSubjElecs in
                %eeg_toolbox
%events         events structure (with 'lfpfile', and 'lfpoffset' fields).
                %see ps2_makeEvents and alignTool
%optional:
%durationMS     %how long to clip
%offsetMS       %where to start relative to mstime
%bufferMS       %buffer to add to beginning and end of clips
%filtfreq       %filters out particular frequencies (e.g., [59.5 60.5]);
                %see buttfilt
%filtorder      %e.g, 4 (see buttfult)

%Note:
%resampleFreq (last arg of gete) is set to 1000 so that gete returns eeg in
%ms. REMOVED
                
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

% Update durationMS and offsetMS with buffer
durationMS = durationMS + 2*bufferMS;
offsetMS = offsetMS-bufferMS;

% replace lfpfile and lfpoffset with eegfile and eegoffset
%First, rename events.eegfile to point to ,noreref
for k=1:length(events)
     events(k).eegfile=events(k).lfpfile;
     events(k).eegoffset=events(k).lfpoffset;
end
    
%check whether it is a bipolar montage or not
if size(elecNum,2) == 1 % not bipolar, proceed normally
   [lfp] = an_gete_ms(elecNum, events,durationMS,offsetMS,bufferMS,...
            filtfreq,filttype,filtorder);
    
elseif size(elecNum,2) == 2 % then it is bipolar
    %First, rename events.eegfile to point to ,noreref
    for k=1:length(events)
         events(k).eegfile=regexprep(events(k).eegfile,'eeg.reref','eeg.noreref');
    end
    
    %Second, an_gete from both channels and subtract them from each other 
    [eeg1] = an_gete_ms(elecNum(1), events,durationMS,offsetMS,bufferMS,...
        filtfreq,filttype,filtorder);
    [eeg2] = an_gete_ms(elecNum(2), events,durationMS,offsetMS,bufferMS,...
        filtfreq,filttype,filtorder);
    
    lfp = eeg1 - eeg2;
end
                
               