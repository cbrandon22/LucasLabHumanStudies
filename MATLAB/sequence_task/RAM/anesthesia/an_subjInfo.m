function[subjInfo,dataDir] = an_subjInfo(subj)
%This function contains the dataDirectory and information for each subject. If a subject name is
%provided as the input, it returns info for only that subject. Otherwise,
%it returns the entire structures

dataDir = '/users/tnl/desktop/c/data/eeg';

i = 1;
subjInfo(i).subj     = 'HUP102_e';
subjInfo(i).sess       = {''}; 
subjInfo(i).rawDataPath  = {'2015-06-24_11-19-32'};
subjInfo(i).logpeg = 0;

i = 2;
subjInfo(i).subj     = 'HUP102_i';
subjInfo(i).sess       = {''}; 
subjInfo(i).rawDataPath  = {'2015-07-06_07-23-48'};
subjInfo(i).logpeg = 0;

i = 3;
subjInfo(i).subj     = 'HUP106_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2015-09-09_12-25-20'};
subjInfo(i).logpeg = 0;

i = 4;
subjInfo(i).subj     = 'HUP108_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2015-10-28_13-55-32'};
subjInfo(i).logpeg = 0;

i = 5;
subjInfo(i).subj     = 'HUP108_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2015-11-06_07-10-39'};
subjInfo(i).logpeg = 0;

i = 6;
subjInfo(i).subj     = 'HUP111_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2016-01-27_20-09-17'};
subjInfo(i).logpeg = 1;

i = 7;
subjInfo(i).subj     = 'HUP114_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2016-03-23_15-56-05'};
subjInfo(i).logpeg = 1;

i = 8;
subjInfo(i).subj     = 'HUP060_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2016-05-06_07-19-48'};
subjInfo(i).logpeg = 1;

i = 9;
subjInfo(i).subj     = 'HUP117_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2016-06-08_17-02-24'};
subjInfo(i).logpeg = 1;

i = 10;
subjInfo(i).subj     = 'HUP116_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2016-05-11_17-16-26'};
subjInfo(i).logpeg = 1;

i = 11;
subjInfo(i).subj     = 'HUP119_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2016-06-24_09-56-00'};
subjInfo(i).logpeg = 1;

i = 12;
subjInfo(i).subj     = 'HUP121_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2016-06-22_13-35-07'};
subjInfo(i).logpeg = 1;

i = 13;
subjInfo(i).subj     = 'HUP121_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2016-07-01_08-19-29'};
subjInfo(i).logpeg = 1;

i = 14;
subjInfo(i).subj     = 'HUP123_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2016-08-15_07-47-01'};
subjInfo(i).logpeg = 1;

% parse input
if ~exist('subj','var') || isempty(subj)
    return
else 
    ind2keep =strcmp({subjInfo.subj},subj);
    if sum(ind2keep) == 0
        error('No Matching Subject Found')
    end
    subjInfo = subjInfo(ind2keep);
end
