function[subjInfo,dataDir] = an_subjInfo(subj)
%This function contains the dataDirectory and information for each subject. If a subject name is
%provided as the input, it returns info for only that subject. Otherwise,
%it returns the entire structures

dataDir = '/users/tnl/desktop/c/data/';

i = 1;
subjInfo(i).subj     = 'HUP119_i';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2016-06-24_09-56-00'};
subjInfo(i).logpeg = 1;

i = 2;
subjInfo(i).subj     = '20150624';
subjInfo(i).sess       = {'pre','post'}; 
subjInfo(i).rawDataPath  = {'2015-06-24_11-19-32'};
subjInfo(i).logpeg = 0;

i = 3;
subjInfo(i).subj     = '20151106';
subjInfo(i).sess       = {'0'}; 
subjInfo(i).rawDataPath  = {'2015-11-06_07-10-39'};
subjInfo(i).logpeg = 0;

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
