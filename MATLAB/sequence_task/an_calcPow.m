function[subjPow] = an_calcPow(subj,elecNum,configNum)
% load bipolar pair of electrodes for entire session
% load events to calculate consciousness points to plot on x axis
subj = 'HUP119_i';
task = 'Consciousness';
evLabel = 'events';
dirs = an_dirs;
dataDir = dirs.data;

% Load config
if ~exist('configNum','var') || isempty(configNum)
    configNum = 1;
end
config = an_config_calcPow(configNum,task);

powDir = fullfile(dirs.scratch,'POWER',num2str(config.configID),subj,[task '_' evLabel]);
cd_mkdir(powDir);
subjPow = [];

%%
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

    [sampleFreq,nBytes,dataformat,gain] = GetRateAndFormat( fileSet{s} );
    sampleFreq = round(sampleFreq);
    
    if isempty(sampleFreq)
        % cannot find params.txt file, cannot determine sampling rate, so error
        error('EEGTOOLBOX:GETE_MS:NOSAMPLERATE','The sample rate for file %s could not be found',fileSet{s});
    end
    
    
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
end

function [samplerate,nBytes,dataformat,gain] = GetRateAndFormat(event)
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
end