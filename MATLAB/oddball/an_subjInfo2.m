function[subjInfo,dataDir] = an_subjInfo2(subj)
%This function contains the dataDirectory and information for each subject. If a subject name is
%provided as the input, it returns info for only that subject. Otherwise,
%it returns the entire structures

dataDir = '/Volumes/HumanStudies/HumanStudies/oddball/eeg';

i = 1;
subjInfo(i).subj     = 'HUP142_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-07-12_12-12-37'};
subjInfo(i).logpeg = 1;

i = 2;
subjInfo(i).subj     = 'HUP143_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-07-17_16-35-55'};
subjInfo(i).logpeg = 1;

i = 3;
subjInfo(i).subj     = 'HUP144_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-07-26_10-37-08'};
subjInfo(i).logpeg = 1;

i = 4;
subjInfo(i).subj     = 'HUP145_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-08-09_11-32-07'};
subjInfo(i).logpeg = 1;

i = 5;
subjInfo(i).subj     = 'HUP147_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-08-28_13-39-47'};
subjInfo(i).logpeg = 1;

i = 6;
subjInfo(i).subj     = 'HUP148_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-08-30_14-21-50'};
subjInfo(i).logpeg = 1;

i = 7;
subjInfo(i).subj     = 'HUP149_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-09-11_11-24-18'};
subjInfo(i).logpeg = 1;

i = 8;
subjInfo(i).subj     = 'HUP150_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-10-03_10-08-19'};
subjInfo(i).logpeg = 1;

i = 9;
subjInfo(i).subj     = 'HUP151_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-10-25_10-18-36'};
subjInfo(i).logpeg = 1;

i = 10;
subjInfo(i).subj     = 'HUP152_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-10-30_10-54-13'};
subjInfo(i).logpeg = 1;

i = 11;
subjInfo(i).subj     = 'HUP153_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-11-22_10-31-08'};
subjInfo(i).logpeg = 1;

i = 12;
subjInfo(i).subj     = 'HUP154_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-11-27_12-28-00'};
subjInfo(i).logpeg = 1;

i = 13;
subjInfo(i).subj     = 'HUP155_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-12-15_19-39-14'};
subjInfo(i).logpeg = 1;

i = 14;
subjInfo(i).subj     = 'HUP156_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-12-22_12-40-11'};
subjInfo(i).logpeg = 1;

i = 15;
subjInfo(i).subj     = 'HUP157_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2017-12-18_12-03-15'};
subjInfo(i).logpeg = 1;

i = 16;
subjInfo(i).subj     = 'HUP159_e';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2018-01-10_10-50-30'};
subjInfo(i).logpeg = 1;

i = 17;
subjInfo(i).subj     = 'HUP166_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2018-03-28_14-18-32'};
subjInfo(i).logpeg = 1;

i = 18;
subjInfo(i).subj     = 'HUP165_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2018-04-03_15-50-12'};
subjInfo(i).logpeg = 1;

i = 19;
subjInfo(i).subj     = 'HUP169_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2018-05-11_07-21-24'};
subjInfo(i).logpeg = 1;

i = 18;
subjInfo(i).subj     = 'HUP171_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2018-06-26_07-19-07'};
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