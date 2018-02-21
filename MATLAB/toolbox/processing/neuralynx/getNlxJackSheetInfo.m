function param_out = getNlxJackSheetInfo(lfpDir,fileExt,parameter,elec)
%
% FUNCTION:
%  getNlxJackSheetInfo.m
%
% DESCRIPTION:
%
% INPUTS:
%  jackSheetFile = ex: 'TJ030';
%  parameter     = path from '/data/eeg/[subj]/raw' to data    
%
% NOTES:
%  (1) written by jfburke 06/12 (john.fred.burke@gmail.com)
% 
% To-Do:
%



% open the jacksheet
jackSheetFileName = sprintf('%s.JacksheetAndParams.txt',fileExt);
jackSheetFile     = fullfile(lfpDir,jackSheetFileName);
if ~exist(jackSheetFile,'file');
  fprintf('\n\n\t!!!EXITING!!!: %s doesn''t exist in %s\n\n',jackSheetFileName,lfpDir)
  return
end
fid=fopen(jackSheetFile,'r');
[X] = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');

% separate teh parmeter names from the values
numParam = length(X);
numElec  = size(X{1},1)-1;
NAMES  = cell(1,numParam);
PARAMS = cell(numElec,numParam);
for k=1:numParam
  NAMES(k)  = X{k}(1);
  PARAMS(:,k) = X{k}(2:end);
end

% get the index of the parameter you want
paramInd = find(strcmp(upper(NAMES),upper(parameter)));

% get the electerode CSC names
elecCSCnames = PARAMS(:,2);

% get the values ofr that parameter
paramVals = PARAMS(:,paramInd);

% get the electrode index
if ~strcmp(upper(elec),'ALL')
  elecInd   = find(strcmp(upper(elecCSCnames),upper(elec)));
  param_all = str2double(paramVals{elecInd});
else    
  param_tmp = nan(numElec,1);
  for k=1:numElec
    param_tmp(k) = str2double(paramVals{k});
  end  
  param_all = unique(param_tmp);
  if length(param_all)>1
    error('more than one param found')
  end
end


param_out=param_all;