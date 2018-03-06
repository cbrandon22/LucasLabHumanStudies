function EEG = geteeg_ms(channel,events,DurationMS,OffsetMS,BufferMS,resampleFreq,RelativeMS)

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
if ~exist('resampleFreq','var')
    resampleFreq = [];
end
if ~exist('RelativeMS','var')
    RelativeMS = [];
end

% Update durationMS and offsetMS with buffer
durationMS = durationMS + 2*bufferMS;
offsetMS = offsetMS-bufferMS;

% reolace lfpfile and lfpoffset with eegfile and eegoffset
%First, rename events.eegfile to point to ,noreref
for k=1:length(events)
     events(k).eegfile=events(k).lfpfile;
     events(k).eegoffset=events(k).lfpoffset;
end

IS_MEF = false;

fileList = {events.eegfile};
fileList(~cellfun( @ischar, fileList )) = {''}; % set any [] filenames to '' for unique to work
fileSet = unique(fileList);

if length(fileSet) <= 1
    if isempty(fileSet) || isempty(fileSet{1})
        EEG = nan(length(events), 0);
        return;
    end
elseif isempty(fileSet{1})
    fileSet = fileSet(2:end);
end
% store the resample filter to save computation
resampleFilt = [];
resampleFiltSamplerate = -1;
resampleExcess = -1;    % extra samples created during resampling

channelSuffix = sprintf('.%03i',channel);
      
if isempty(resampleFreq)
    % if resampleFreq is not specified, always re-sample to the first
    % session's sampling rate to ensure the same number of samples for all
    % events, across sessions
    resampleFreq = round(GetRateAndFormat2( fileSet{1} ));
else
    resampleFreq = round(resampleFreq);
end

% base final datasize on resampled data
buffDur = ceil( (DurationMS+2*BufferMS)*resampleFreq/1000 );
buffer = fix( (BufferMS)*resampleFreq/1000 );
duration = buffDur - 2*buffer;

% allocate space for EEG data
EEG = nan(length(events),duration);

% loop through the files
for s = 1:length(fileSet)
    fileMask = strcmp(fileSet{s}, fileList);
    
    [sampleFreq,nBytes,dataformat,gain] = GetRateAndFormat2( fileSet{s} );
%     sampleFreq = eegparams(samplerate,fileSet{s});
    sampleFreq = round(sampleFreq);
    
    
    if isempty(sampleFreq)
        % cannot find params.txt file, cannot determine sampling rate, so error
        error('EEGTOOLBOX:GETE_MS:NOSAMPLERATE','The sample rate for file %s could not be found',fileSet{s});
    end
    
        %compute input params for gete based on next session's sampling rate
    readDuration = ceil( (DurationMS+2*BufferMS)*sampleFreq/1000 );
    readOffset = round( (OffsetMS-BufferMS)*sampleFreq/1000 );
    
    %%% open file and read in EEG
    eegfname = [fileSet{s}, channelSuffix];
    startIndices = [events(fileMask).eegoffset] + readOffset;
    
    fileEEG = nan(length(startIndices),max(readDuration,buffDur));
    
    if IS_MEF
        [readEEG,fileEEGlen] = decomp_mef_events(eegfname,startIndices,readDuration,'');
        fileEEG(:,1:readDuration) = readEEG';
        % EH - need to add code (either here or to decomp_mef_events) to
        % check the 'unpadded lead' filename (see a few lines below)
    else
        fileEEGlen = zeros(length(startIndices),1,'uint32');
        
        eegfile = fopen(eegfname,'r','l'); % NOTE: the 'l' means that it came from a PC!
        
        if eegfile == -1
            eegfname = sprintf('%s.%i',fileSet{s},channel); % now try unpadded lead
            eegfile = fopen(eegfname,'r','l');
        end

        % if the file still can't be opened, throw an error
        if eegfile == -1                
            error('Missing EEG file: %s\n',eegfname);
        end
        
        for i = 1:length(startIndices)
            status = fseek(eegfile,nBytes*startIndices(i),-1);    
            if status == 0
                readbytes = fread(eegfile,readDuration,dataformat)';
                fileEEG(i,1:length(readbytes)) = double(readbytes);
                fileEEGlen(i) = length(readbytes);
            end
        end
        
        fclose(eegfile);
    end
    
    noEEGmask = (fileEEGlen == 0);
    fullEEGmask = (fileEEGlen == readDuration);
    shortEEGind = find( ~(fullEEGmask | noEEGmask) );
    
    if any(noEEGmask)
        warning('EEGTOOLBOX:GETE_MS:NODATA','%s: EEG data were not found for %d event(s) -- setting NaN',...
                eegfname,sum(noEEGmask));
    end
    
    if ~isempty(shortEEGind)
        warning('EEGTOOLBOX:GETE_MS:INCOMPLETEDATA','%s: not all samples read for %d event(s) -- appending NaN',...
                eegfname,length(shortEEGind));
    end
    
    if IS_MEF
        % set unread EEG data to NaN
        for i = 1:length(shortEEGind)
            ind = shortEEGind(i);
            fileEEG(ind,fileEEGlen(ind)+1:end) = NaN;
        end
        fileEEG(noEEGmask,:) = NaN;
    end
    
    if( resampleFreq ~= sampleFreq )
        if( resampleFiltSamplerate ~= sampleFreq )
            [resamplebytes,resampleFilt] = resample(fileEEG(fullEEGmask,1:readDuration)',resampleFreq,sampleFreq);
            resampleFiltSamplerate = sampleFreq;
            resampleExcess = fix( (size(resamplebytes,1) - buffDur)/2 );
        else
            resamplebytes = resample(fileEEG(fullEEGmask,1:readDuration)',resampleFreq,sampleFreq,resampleFilt);
        end
        
        % often, resamplebytes will be longer than buffDur
        % b/c resample uses: length(output) = ceil(length(input)*resampleFreq/sampleFreq)
        fileEEG(fullEEGmask,1:buffDur) = resamplebytes(resampleExcess+1:resampleExcess+buffDur,:)';

        for i = 1:length(shortEEGind)
            ind = shortEEGind(i);
            % no need to check that resampleFilt is the correct sampling rate, since that was just done above
            resamplebytes = resample(fileEEG(ind,1:fileEEGlen(ind))',resampleFreq,sampleFreq,resampleFilt);
            fileEEG(ind,1:length(resamplebytes)) = resamplebytes';
        end
    end
    
    EEG(fileMask,:) = gain .* fileEEG(:,buffer+1:duration+buffer);
end

%JRM NOTE: this is old code -- I left it in from the old version of gete_ms
%take relative baseline correction
if length(RelativeMS) == 2
    % get the average for the range
    offset = fix((OffsetMS)*resampleFreq/1000);
    relative = fix((RelativeMS)*resampleFreq/1000);
    relative = relative - offset + 1;
    relative(2) = relative(2) - 1;

    % calculate the relative
    releeg = mean(EEG(:,relative(1):relative(2)),2);

    % subtract the baseline
    EEG = bsxfun(@minus, EEG, releeg);
end

function [samplerate,nBytes,dataformat,gain] = GetRateAndFormat2(event)
%GETRATEANDFORMAT - Get the samplerate, gain, and format of eeg data.
%
% function [samplerate,nBytes,dataformat,gain] = GetRateAndFormat(event)
%

if ischar(event)
    % event is actually a path to an EEG file
    path = event;
else
    path = event.eegfile;
end

% Look in "noreref" directory for session-specific parameters file
fs = filesep;
sepInds = strfind(path,fs);

if isempty(sepInds) % relative path (current folder)
    dir = '';
    parentDir = ['..' fs]; % parent dir is up one level
    name = path;
elseif length(sepInds) == 1 % relative path (child folder)
    dir = path(1:sepInds(end));
    parentDir = ''; % parent dir is current dir
    name = path(sepInds(end)+1:end);
else
    dir = path(1:sepInds(end));
    parentDir = path(1:sepInds(end-1));
    name = path(sepInds(end)+1:end);
end

nameDotInds = strfind(name,'.');
if ~isempty(nameDotInds) && nameDotInds(1) ~= 1
    baseName = name(1:nameDotInds(1)-1); % strip all extensions
else
    baseName = name;
end

sessionParamsPath = [parentDir 'lfp.noreref' fs 'params.txt']; % removed basename from path (e.g. HUP119_i_24Jun16_0959)

% EH - trying to open the file is a quicker way to check for its existence
% than using exist(sessionParamsPath,'file')
file = fopen(sessionParamsPath,'rt'); 
if( file == -1 )
    sessionParamsPath = [dir 'params.txt'];
    file = fopen(sessionParamsPath,'rt'); 
end

params = eegparams({'samplerate','gain','dataformat'},file);

samplerate = params{1};
if( isempty(samplerate) )
    % EH - can't do anything if the sample rate isn't present (no default)
    gain = [];
    dataformat = '';
    nBytes = [];
    return;
end

if( ~isempty(params{2}) )
    gain = params{2};
else
    gain = 1.0;
end

if( ~isempty(params{3}) )
    dataformat = params{3};
else
    dataformat = 'short';
end

switch dataformat
    case {'short','int16'}
        nBytes = 2;
    case {'single','int32'}
        nBytes = 4;
    case {'double','int64'}
        nBytes = 8;
    otherwise
        error('BAD DATA FORMAT!');
end

function params = eegparams(field,filepath)
%EEGPARAMS - Get a subject specific eeg parameter from the params.txt file.
% 
% If paramdir is not specified, the function looks in the 'docs/'
% directory for the params.txt file.
%
% The params.txt file can contain many types of parameters and will
% evaluate them as one per line.  These are examples:
%
% Channels 1:64
% samplerate 256
% subj 'BR015'
%
% FUNCTION:
%   p = eegparams(field,paramdir)
%
% INPUT ARGS:
%   field = 'samplerate';        % Field or cell array of fields to retrieve
%   filepath = '~/eeg/012/dat/params.txt';  % Path to parameter file.
%       Alternatively, filepath can be a handle to an open file
%
% OUTPUT ARGS:
%   params- the parameters in a cell array, evaluated with eval()
%

if( ~iscell(field) )
    noCell = true;
    field = {field};
else
    noCell = false;
end

if( ~ischar(filepath) )
    file = filepath; % filepath is actually a file handle
else
    file = fopen(filepath,'rt'); % open the file
end

if( file == -1 )
    fileC = {'',''};
else
    fileC = textscan(file,'%s%s','Delimiter','= \b\t');
    fclose(file);
end

nField = length(field);
params = cell(nField,1);

for i = 1:nField
    ind = find( strcmp(field{i},fileC{1}) );
    
    if( ~isempty( ind ) )
        params{i} = eval(fileC{2}{ind(end)});
    else
        params{i} = [];
    end
end

if( noCell )
    % since input was not cell array, output shouldn't be either
    params = params{1};
end