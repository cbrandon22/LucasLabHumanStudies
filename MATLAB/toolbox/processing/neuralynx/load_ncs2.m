function [data,info,TSblock,TSsample]=load_ncs2(filename)
%
% FUNCTION:
%   [data,info,allTS]=load_ncs2(filename)
%
% DESCRIPTION:
%   this function returns the raw 2 byte samples inside this .ncs
%    
% INPUTS:
%   filename: path to the ncs file
%
% OUTPUTS:
%   data:  LFP data in microvolts 
%   info:  the output of stat_ncs.m (eeg_toolbox under basic)
%   TSblock: timestamps of each block 
%   TSsample: timestamps of each sample (interpolated)
%
% NOTES:
%  (1) written 08/11
%  (2) It is the second version of load_ncs.m. Differences incude:
%       - that it returns a timestamp for each sample, not just for
%         each 512 sample block of data.   
%       - automatically multiplies the data by the gain
%
% To-Do:
%   (1) Should we move this to the eeg_toolbox?
%   (2) Move stat_ncs to the units toolbox?
%   (3) Shold we write out a downsampled version of the data and
%       save it in eeg.noreef?
%

samplesPerRecord=512;
bytesPerRecord=1024+20;
headerBytes=16*1024;

info=stat_ncs(filename);
disp(sprintf('difference between samples is %g microseconds',info.sampleTSDiff));
disp(sprintf('header sampling rate is %g Hz',info.headerSampleRate));
disp(sprintf('actual sampling rate is %g Hz',info.actualSampleRate));

if info.numPotentialSamples~=info.numActualSamples
  warning('missing some samples...');
end

% preallocation
data     = zeros(1,info.numPotentialSamples,'int16');
TSblock  = zeros(1,info.numRecords);
TSsample = zeros(1,info.numPotentialSamples);
prevTS   = nan;

% open the file
fid=fopen(filename,'r','l');
fseek(fid,headerBytes,'bof'); 

% loop thru all blocks of 512 samples
ticker=0;
tick_inc=10;
fprintf('Loading data:')
for record=0:info.numRecords-1
  if record./(info.numRecords-1)*100 >= ticker
    fprintf(' %d%%',ticker)
    ticker=ticker+tick_inc;
  end
  
  % the timestamp this block
  ts=fread(fid,1,'int64');

  % Three other other things for this blocks
  %  1. 
  %  2. 
  %  3. number valid sample
  % read in 3 in one shot for speed
  tmp=fread(fid,3,'int32');
  numValidSamp=tmp(3);
  
  thisRecordData=fread(fid,512,'*int16');
  firstIndex=round((ts-info.firstSampleTime)/info.sampleTSDiff);
  i=1:numValidSamp;
  data(firstIndex+i)=thisRecordData(i);

  if ~isnan(prevTS)& abs(ts-(prevTS+info.sampleTSDiff*512))>10 
    %if difference is greater than 10 microseconds
    %note: 10 microseconds is totally arbitrary... someone should
    %think about this more.  i guess this shit is due to rounding,
    %as indicated by some .pdf in the neuralynx documentation, but
    %this does seem like a rather large difference....   
    warning('timestamp difference has changed... yuck');
  elseif (numValidSamp~=samplesPerRecord) %record < info.numRecords-1 &&
      if record < info.numRecords-1
          warning('512 samples WERE NOT read. WTF');
      else
          continue; % don't display warning if last block
      end
  else
  prevTS=ts;
  TSblock(record+1)=ts;
  TSsample(firstIndex+i)=ts+[0:info.sampleTSDiff:info.sampleTSDiff*(samplesPerRecord-1)];
  end
  % perform this check just to be sure
  %if ~isequal(TSsample(1:samplesPerRecord:samplesPerRecord*record+1),TSblock(1:record+1))
  %  error('something went bad with your interpolation')
  %end
  
  % another check
  %if sum(unique(diff(TSsample(1:512*record+1)))-info.sampleTSDiff>10)~=0;
  %  error('something went bad with your interpolation')
  %end
end
data=data.*(info.bitVolts*1e6);
fprintf('\n\n')
fclose(fid);