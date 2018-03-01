function getNlxHeaderFileParams(chanFile,bytesInHeader)
%
% FUNCTION:
%  getNlxHeaderFileParams.m
%
% DESCRIPTION:
%
% INPUTS:
%
% NOTES:
%  (1) written by jfburke 05/12 (john.fred.burke@gmail.com)
%

if ~exist('bytesInHeader','var')||isempty(bytesInHeader)
  bytesInHeader=16*1024;
end

%read the header into a string
fid=fopen(chanFile,'r','l');

% separate the header into the commented and actual header
commentedHeader = [];
actualHeader    = [];
while true
  N=fgetl(fid);
  if N==-1;break;end;
  if strcmp(N(1),'#')
    commentedHeader(end+1,:)=N;
  else
    actualHeader(end+1,:)=N;
  end
end
keyboard