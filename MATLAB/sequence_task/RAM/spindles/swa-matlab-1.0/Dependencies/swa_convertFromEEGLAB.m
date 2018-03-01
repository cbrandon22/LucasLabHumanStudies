function [Data, Info] = swa_convertFromEEGLAB(fileName, filePath)

% ask the user for the eeg .set file
if nargin < 1
    [fileName, filePath] = uigetfile('*.set');
end

% user output
fprintf(1, 'loading set file: %s...', fileName);
load(fullfile(filePath, fileName), '-mat');
fprintf(1, 'done \n');

% TODO: make loading compatible with old .dat files from EEGLAB
if strcmp(EEG.data(end-3: end), '.dat')
    fprintf(1, 'Warning: currently not compatible with .dat files, use EEGLAB to convert to .fdt files');
    return;
end

% allocate the data information to the Info structure
Info.Recording.dataFile   = EEG.data;
Info.Recording.dataDim    =[EEG.nbchan, EEG.pnts];
Info.Recording.sRate      = EEG.srate;
Info.Recording.reference  = EEG.ref;
Info.Electrodes           = EEG.chanlocs;

% open the binary file and read contents to the Data.Raw structure
fprintf(1, 'loading data file: %s...', fileName);
fid = fopen(fullfile(filePath, Info.Recording.dataFile));
Data.Raw = fread(fid, Info.Recording.dataDim, 'single');
fclose(fid);
fprintf(1, 'done \n');
