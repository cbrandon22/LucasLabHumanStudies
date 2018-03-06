function pyFR_RECOG_raster(subj,rawDataPath,behDataPath,channelStr)
%
% FUNCTION:
%   
% DESCRIPTION:
%
% INPUT:
%
% OUTPUT:
%
% NOTES:
%   (1) written by jfburke [''] (john.fred.burke@gmail.com)
%
%

clear global
global rootDir CHANNELSTR rawDir epoch win

CHANNELSTR = channelStr;
if ~ismac
  rootDir = '';
else
  rootDir = '/Volumes/RHINO_root';
end

win     = [250];

dataDir = fullfile(rootDir,'/data/eeg/');
subjDir = fullfile(dataDir,subj);
rawDir  = fullfile(subjDir,'raw',rawDataPath);
behDir  = fullfile(subjDir,'behavioral/pyFR',behDataPath);

% load events
evFile  = fullfile(behDir,'RECOG_events.mat');
if ~exist(evFile,'file')
  error('no events data')
end
ev = load(evFile);
ev = ev.events;

% separate out recalled and non-recalled events
respEv = ev(strcmp({ev.type},'RECOG_RESP'));
presEv = ev(strcmp({ev.type},'RECOG_PRES'));

epoch   = [-1000 5000];
oldEv = presEv([presEv.isTarget]==1);
newEv = presEv([presEv.isTarget]==0);
mkRaster_local(oldEv,newEv,1)

epoch  = [-5000 2000];
tar_Ev = respEv([respEv.isTarget]==1);
lur_Ev = respEv([respEv.isTarget]==0);
mkRaster_local(tar_Ev,lur_Ev,4)

yes_Ev = respEv([respEv.recognized]==1);
no_Ev  = respEv([respEv.recognized]==0);
mkRaster_local(yes_Ev,no_Ev,7)

hit_Ev  = respEv([respEv.isTarget]==1 & [respEv.recognized]==1);
miss_Ev = respEv([respEv.isTarget]==1 & [respEv.recognized]==0);
mkRaster_local(hit_Ev,miss_Ev,10)

fa_Ev  = respEv([respEv.isTarget]==0 & [respEv.recognized]==1);
cr_Ev = respEv([respEv.isTarget]==0 & [respEv.recognized]==0);
mkRaster_local(fa_Ev,cr_Ev,13)


% get the spieks times
function mkRaster_local(recEv,notEv,fNum)

  global rootDir CHANNELSTR rawDir epoch win

  % get lfp file info
  [lfpDir,fileExt] = fileparts(recEv(1).lfpfile);
  lfpDir           = fullfile(rootDir,lfpDir);

  timesFileName = sprintf('times_%s.mat',CHANNELSTR);
  timesFile     = fullfile(rawDir,timesFileName);
  if ~exist(timesFile,'file')
    error('no times data')
  end
  tData    = load(timesFile);
  cc       = tData.cluster_class;
  TS       = getNlxJackSheetInfo(lfpDir,fileExt,'DOWNSAMPLE_TSDIFF',[CHANNELSTR '.ncs']);
  
  % get unique number of clusters
  unClusters_all = unique(cc(:,1));
  unClusters     = unClusters_all(unClusters_all~=0);
  if isempty(unClusters)
    fprintf('no clusters\n\n')
    return
  end
  
  % loop through clusters
  for k=1:length(unClusters)
    thisC    = unClusters(k);
    spkTimes = cc([cc(:,1)]==thisC,2);
    
    % get the raster
    [rastRec,tVect] = getRastAndFr_local(spkTimes,[recEv.lfpmstime],epoch,TS);
    [rastNot,tVect] = getRastAndFr_local(spkTimes,[notEv.lfpmstime],epoch,TS);
    
    % plot the rasters
    axRec = plotRaster(rastRec,tVect,'r',TS,1:length(recEv),fNum);
    axNot = plotRaster(rastNot,tVect,'b',TS,1:length(notEv),fNum+1);
    
    % get the smooth firing rate
    numSamInWin = floor(win./(TS*1000));
    MA          = 1/numSamInWin*[1:numSamInWin];
    frRec       = filter(MA,1,rastRec,[],2)*(1000./win);
    frNot       = filter(MA,1,rastNot,[],2)*(1000./win);
    
    % plot the  
    figure(fNum+2);clf  
    ci_tmp = std(frNot,[],1)./sqrt(size(frNot,1));
    errorbar(tVect,mean(frNot,1),ci_tmp,'b');    
    hold on  
    ci_tmp = std(frRec,[],1)./sqrt(size(frRec,1));
    errorbar(tVect,mean(frRec,1),ci_tmp,'r');  
    yL=get(gca,'YLim');
    plot(zeros(1,100),linspace(yL(1),yL(2),100),'--k','LineWidth',2)
    formatPic(gca)
    hold off
    xlim([tVect(1)+win tVect(end)])
    
    % test the vetors  
    fprintf('Cluster %d. pausing.....',k)
    pause
    fprintf('done\n')
  end
  
function ax=plotRaster(rast,tVect,col,Ts,sp,fNum)
  [y,x]=find(rast);
  TsMSec   = Ts*1000;   
  xT       = 500;
  TsSam    = floor(xT./TsMSec);
  figure(fNum)
  clf
  ax=axes;
  cla(ax);
  set(ax,'XTick',1:TsSam:size(rast,2),'YTick',1:4:size(rast,1));  
  set(ax,'XTIckLabel',tVect(1):xT:tVect(end))
  set(ax,'YTIckLabel',1:4:length(sp))  
  set(ax,'YLim',[0 size(rast,1)+1])
  line('Parent',ax,'XData',x,'YData',y,'Marker','o',...
       'MarkerFaceColor',col,'MarkerEdgeColor','none','MarkerSize',4,...
       'LineStyle','none');  
  xlabel('Times (ms)')
  ylabel('Trial')
  %formatPic(ax)
  
function [rast,timVect]=getRastAndFr_local(spTimes,evTimes,epoch,Ts);
  TsMSec      = Ts*1000;   
  numSamInTim = floor(diff(epoch)./TsMSec);    
  timVect     = linspace(epoch(1),epoch(2),numSamInTim);

  % allocate 
  nEv   = length(evTimes);
  rast  = false(nEv,numSamInTim);
  
  % loop thru events
  for k=1:nEv
    thisTimBound    = epoch+evTimes(k);    
    spTimThisEv_abs = spTimes(spTimes>=thisTimBound(1) & spTimes<=thisTimBound(end));
    spTimThisEv     = spTimThisEv_abs-evTimes(k);
    if isempty(spTimThisEv)
      continue
    end
    nSp = length(spTimThisEv);
    
    % find closest index to actual time
    M            = repmat(timVect',1,nSp)-repmat(spTimThisEv',numSamInTim,1);
    [~,minInd]   = min(abs(M),[],1);
    rast(k,minInd) = true;
  end

  
  
