function SpikeLFP_viewer(subj,rawDataPath,channelStr,labelSpikes,notchIt)
%
% FUNCTION:
%   plotSpikeAndLFPTimeSeries(subj,rawDataPath,channelStr,doFilt)
%
% DESCRIPTION:
%   plots time series of spikes with LFP.  Slides window with 90%
%   overalp.  This function is a stripped down version of
%   spikeBrowser.m written by Joshua Jacobs.  Wraps around
%   load_ncs.m 
%
% INPUTS:
%  subj= ex: 'TJ030';
%  rawDataPath= path from '/data/eeg/[subj]/raw' to data    
%  channelStr= ex: 'CSC129';  
%  doFilt=false
%
% NOTES:
%  (1) written by jfburke 08/11 (john.fred.burke@gmail.com)
%  (2) assumes ncs files are of form [channelStr].ncs
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADDITIONAL CONFIGURABLE PARAMS
OVERLAP  = .9; % percatage to overlap windows
clusters = 'ALL'; % clusters to plot  
winLen   = 5000; % size of window in ms
incLen   = floor(winLen*OVERLAP);
msecBetweenSpikes = 3;
msecAroundSpikes  = 3;
YL      = [-100 100];
YL_raw  = [-600 600];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure out what type of data we're using
recType = 'cheetah';
if strfind(upper(channelStr),'NSX')
  thisChanFile = sprintf('%s.NC5',channelStr);    
  recType = 'black';
elseif ~strcmp(upper(channelStr),'EVENTS')
  thisChanFile = sprintf('%s.ncs',channelStr);
else strcmp(upper(channelStr),'EVENTS')
  thisChanFile = sprintf('%s.nev',lower(channelStr));
end

% make the directories where the data live
if ~ismac
    if strcmp(recType,'cheetah')
        dataDir    = '/data/eeg';
        subjDir    = fullfile(dataDir,subj);
        rawDir     = fullfile(subjDir,'raw');
        rawDataDir = fullfile(rawDir,rawDataPath);
    else
        dataDir    = '/data/continuous';
        subjDir    = fullfile(dataDir,subj);
        rawDataDir = fullfile(subjDir,rawDataPath);
    end
else
    dataDir    = '/Volumes/RHINO_root/data/eeg';
    subjDir    = fullfile(dataDir,subj);
    rawDir     = fullfile(subjDir,'raw');
    rawDataDir = fullfile(rawDir,rawDataPath);
end
timesFile     = sprintf('times_%s.mat',channelStr);
thisChanPath  = fullfile(rawDataDir,thisChanFile);
timesFilePath = fullfile(rawDataDir,timesFile);

% load the data
if strcmp(recType,'cheetah')
    [data,info,TSblock,TSsamples]=load_ncs2(thisChanPath);
    sr = info.actualSampleRate;
    sampleTSDiff = info.sampleTSDiff;
    samPerInc  = ceil(incLen/sampleTSDiff*1000);
    samPerWin  = ceil(winLen/sampleTSDiff*1000);
else
    [data,info] = load_nc5(thisChanPath);    
    sr = info.sr;
    TSsamples = info.TimeStamps;
    samPerInc = ceil(incLen/1000*sr);
    samPerWin = ceil(winLen/1000*sr);
end

% get looping params
incPerFile = floor((length(data)-samPerWin)./samPerInc);

% set up the plotting tools
close all;
[ax ax_raw LFP LFP_raw zeroLine zeroLine_raw] = mkThisFig_data(1,YL,YL_raw,samPerWin);
%[ax_spk] = mkThisFig(2,YL,samPerWin);ax_spk=[];
%[ax_raw] = mkThisFig(3,YL,samPerWin);
ax_spk=[];

% get the times of all spikesa
sTimes=getSpikeTimes_local(rawDataDir,channelStr); 
LoopTheseClusters=1:length(sTimes);

% directions:
fprintf('\n\nPress ''enter'' to go forward.\n')
fprintf('Press ''b'' and then ''enter'' to go backward.\n')
fprintf('Looping....\n')

% loop thru and plot
spikes = [];
count = 0;
for k=1:incPerFile
  if count>incPerFile;break;end
  count=count+1;
  thisWindowSamStart=(count-1)*samPerInc;
  slidingWindowInd = thisWindowSamStart+1:thisWindowSamStart+samPerWin; 
  thisTS  = TSsamples(slidingWindowInd)/1000;% - TSsamples(1)/1000; % now in msec
  thisLFP_raw = double(data(slidingWindowInd));  
  thisLFP = buttfilt(thisLFP_raw,[500 9000],sr,'bandpass',4);
  %thisLFP = buttfilt(thisLFP_raw,[2 10],round(info.actualSampleRate),'bandpass',2);
  if notchIt
    thisLFP_raw = buttfilt(thisLFP_raw,[58 62],round(sr),'stop',2); 
    thisLFP_raw = buttfilt(thisLFP_raw,[2 300],round(sr),'bandpass',2); 
    %thisLFP_raw = buttfilt(thisLFP_raw,[2 10],round(info.actualSampleRate),'bandpass',2);
  end    
  
  % plot the LFP signal
  set(LFP,'YData',thisLFP,'XData',thisTS-thisTS(1));  
  set(LFP_raw,'YData',thisLFP_raw,'XData',thisTS-thisTS(1));  
  
  if sum(isnan(thisTS))>0
    set(ax,'xlim',[thisTS(1)-thisTS(1) thisTS(end)-thisTS(1)])
    %set(ax_raw,'xlim',[thisTS(1) thisTS(end)])
  end
  set(zeroLine,'XData',thisTS-thisTS(1))
  set(zeroLine_raw,'XData',thisTS-thisTS(1))
  
  % get the times of all spikes
  %if labelSpikes
  %  [dots dots_raw] = plotSpikeTimes_local(LoopTheseClusters,sTimes,...
%				  thisTS,ax,ax_raw);			 
  %else
  %  dots=[];
  %end
  
  if mod(k,10)==0
    fprintf('%d%% ',round(k/incPerFile*100))
    if mod(k,100)==0
      fprintf('\n')
    end
  end
  
  in = input('','s');
  if strcmp(upper(in),'B')
    fprintf('going back on window\n')
    count=count-2;
  end
  
  % remove the dots
%   for k=1:length(dots)
%     if k==1
%       pause(.5)
%     end
%     delete(dots(k));
%     delete(dots_raw(k));
%   end
end

size(spikes)

function [l1 l2]= plotSpikeTimes_local(LTC,sT,t,ax,ax_raw); 
  
  spikeCols  = {'b','r','g','c','m','y'};  
  l1 = [];l2=[];
  for cInd=1:length(LTC)
    c=LTC(cInd);
    SpikeIndThisWin=find(sT{c}>t(1) & sT{c}<=t(end));
    for s=1:length(SpikeIndThisWin)
      thisSpike = sT{c}(SpikeIndThisWin(s));
      %line('Parent',ax,'YData',linspace(YL(1),YL(2),100),'XData',ones(1,100)*thisSpike,...
      % 'Color',spikeCols{c},'LineStyle','--','LineWidth',1);
      %line('Parent',ax,'YData',YL(2)-YL(2)*.25,'XData',thisSpike,...
      % 'Color',spikeCols{c},'Marker','.','MarkerSize',30);
      l1_tmp = line('Parent',ax,'YData',-10,'XData',thisSpike-t(1),...
		    'Color',spikeCols{c},'Marker','.','MarkerSize',30);
      l2_tmp = line('Parent',ax_raw,'YData',0,'XData',thisSpike-t(1),...
		    'Color',spikeCols{c},'Marker','.','MarkerSize',30);
      l1 = [l1; l1_tmp];
      l2 = [l2; l2_tmp];
    end    
  end  

function sTimes=getSpikeTimes_local(rawDataDir,channStr); 
thisTimesFile = sprintf('times_%s.mat',channStr);
thisTimesPath = fullfile(rawDataDir,thisTimesFile);

% does it exist?
f=fopen(thisTimesPath);
if f==-1;sTimes=[];return;end;fclose(f);

% get the cluster_class variable
X=load(thisTimesPath);
cc=X.cluster_class;
numSpikes=max(cc(:,1));

% organize spikes, skip cluster 0
sTimes=cell(1,numSpikes);
for k=1:numSpikes
  sTimes{k}=cc(cc(:,1)==k,2);
end


function [l spikes]= getSpikeTimes_local_foo(x,x_raw,t,ax,ax_raw,ax2,timeMax,spikeWin,dt,raw); 
  % get indices of spikes
  spikeInd = find(x>th);    
  if isempty(spikeInd)
    l=[];spikes=[];return;
  end
  N            = length(x);
  spikeSamVect = -spikeWin:spikeWin;
  spikeTimVect = spikeSamVect*dt;
  
  % find the maximum value for each spike (separated by timeMax msec)
  firstSpikes  = [0 find((diff(spikeInd,1)>timeMax)) length(spikeInd)];
  centerSpikes = []; 
  for fs=1:length(firstSpikes)-1
    theseInd         = spikeInd(firstSpikes(fs)+1:firstSpikes(fs+1));
    [~,ind_foo]      = max(x(theseInd));
    centerSpikes_tmp = theseInd(ind_foo);
    if (centerSpikes_tmp >= (N-spikeWin)) | (N <= spikeWin)
      continue
    end
    centerSpikes = cat(1,centerSpikes,centerSpikes_tmp);
  end
  
  % clip spikeWin msec around each spike 
  spikes   = nan(length(centerSpikes),2*spikeWin+1);
  spikeRaw = nan(length(centerSpikes),2*spikeWin+1);
  for cs=1:length(centerSpikes)
    thisCenter = centerSpikes(cs);
    thisSpike  = x(thisCenter-spikeWin:thisCenter+spikeWin);
    spikes(cs,:) = thisSpike;
  %  
  %  % plot it
  %  line('Parent',ax2,'XData',spikeTimVect,'YData',thisSpike,...
  %	 'LineWidth',2);
  %    set(ax2,'XLim',[spikeTimVect(1) spikeTimVect(end)]);
  %    set(ax2,'YLim',3*[-th th]);
    end
    
  % plot the spikes on the axes with the sliding window
  l = nan(1,length(spikeInd));
  for s=1:length(spikeInd)
    l(s)=line('Parent',ax,'YData',th,'XData',t(spikeInd(s)),...
	 'Color','r','Marker','.','MarkerSize',30);
  end    

function [ax1 ax2 L_1 L_2 zL_1 zL_2] = mkThisFig_data(f,YL_1,YL_2,samPerWin)
  figure(f);
  set(f,'units','normalized','Position',[0 0 1 .7]);
  ax1=subplot(2,1,1);
  ax2=subplot(2,1,2);
  set(ax1,'ylim',YL_1,'LineWidth',3,'FontWeight','bold','FontSize',15,'FontName','courier');
  set(ax2,'ylim',YL_2,'LineWidth',3,'FontWeight','bold','FontSize',15,'FontName','courier');

  zL_1=line('Parent',ax1,'YData',zeros(1,samPerWin),'XData',zeros(1,samPerWin),...
	    'LineWidth',2,'LineStyle','--','Color','k');
  zL_2=line('Parent',ax2,'YData',zeros(1,samPerWin),'XData',zeros(1,samPerWin),...
	    'LineWidth',2,'LineStyle','--','Color','k');
  
  L_1=line('Parent',ax1,'YData',zeros(1,samPerWin),'XData',zeros(1,samPerWin),...
	   'LineWidth',1,'LineStyle','-','Color','k');
  L_2=line('Parent',ax2,'YData',zeros(1,samPerWin),'XData',zeros(1,samPerWin),...
	   'LineWidth',1,'LineStyle','-','Color','k');
  
  xlabel(ax1,'Time (ms)');
  ylabel(ax1,'Potential (\muV)');
  
  xlabel(ax2,'Time (ms)');
  ylabel(ax2,'Potential (\muV)');
  
function [ax L zL] = mkThisFig(f,YL,samPerWin)
  figure(f);ax=axes;
  set(ax,'ylim',YL,'LineWidth',3,'FontWeight','bold','FontSize',15,'FontName','courier')
  zL=line('Parent',ax,'YData',zeros(1,samPerWin),'XData',zeros(1,samPerWin),...
	  'LineWidth',2,'LineStyle','--','Color','k');
  L=line('Parent',ax,'YData',zeros(1,samPerWin),'XData',zeros(1,samPerWin),...
	 'LineWidth',1,'LineStyle','-','Color','k')
  xlabel(ax,'Time (ms)')
  ylabel(ax,'Potential (\muV)')
