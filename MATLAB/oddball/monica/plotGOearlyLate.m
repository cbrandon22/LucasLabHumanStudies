config = le_config_calcPow(1,'oddball');
dirs = le_dirs;
subj = 'HUP142_i';
fileExt = 'HUP142_i_12Jul17_1219';
response = 'go';

load('/Users/monicakhattak/Desktop/LucasLAB/data/events/oddball/HUP142_i_events.mat')
% 
% open the jacksheet
jackSheetFileName = sprintf('%s.JacksheetAndParams.txt',fileExt);
jackSheetFile     =  fullfile(dirs.data,'eeg','oddball',subj,'lfp.noreref',jackSheetFileName);
if ~exist(jackSheetFile,'file');
  fprintf('\n\n\t!!!EXITING!!!: %s doesn''t exist in %s\n\n',jackSheetFileName,lfpDir)
  return
end
fid=fopen(jackSheetFile,'r');
[X] = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
%[X] = textscan(fid,'%d%s','delimiter','\t');

% separate teh parmeter names from the values
numParam = length(X);
numElec  = size(X{1},1)-1;
NAMES  = cell(1,numParam);
PARAMS = cell(numElec,numParam);
for k=1:numParam
  NAMES(k)  = X{k}(1);
  PARAMS(:,k) = X{k}(2:end);
end

ElecNum = PARAMS(:,1);
elecs = str2double(ElecNum);
ElecName = PARAMS(:,2);
Enames = char(ElecName);

%[passEvents] = splitevents(events);
[retain1, retain2] = vplotter(events, 'goEARLY');
     
% %OPTIONS:
% %backgroundRESPONSE
% %goNOGO
% %responseERRORS
% %gonogoERRORS
     
    % if an explant:
%      earlyTarg = passEvents(retain1);
%      earlyBack = passEvents(retain2);
     %If an implant:
     earlyTarg = events(retain1);
     earlyBack = events(retain2);
% 
    eTargEEG = zeros(sum(retain1),config.durationMS,length(elecs)); %%%
    eBackEEG = zeros(sum(retain2),config.durationMS,length(elecs));
   [ev_out]=get_baseline_random(events,config.baseSecBetEv,config.baseJitt); %%%
    %%%change rec_type to lfp but after remember to change it to eeg



for i = 1:length(elecs)
    eeg1 = an_getlfp_ms_wrapper (elecs(i),earlyTarg,config.durationMS, config.offsetMS, config.bufferMS, (50.5),'low', 14);
    eeg2 = an_getlfp_ms_wrapper (elecs(i),earlyBack,config.durationMS, config.offsetMS, config.bufferMS, (50.5),'low', 14 );
    eeg3 = an_getlfp_ms_wrapper (elecs(i),ev_out,config.durationMS, config.offsetMS, config.bufferMS, (50.5),'low', 14 );
    baseEEG(:,:,i) = eeg3(:,config.bufferMS+1:(size(eeg3,2)-config.bufferMS));
    eBackEEG(:,:,i) = eeg2(:,config.bufferMS+1:(size(eeg1,2)-config.bufferMS));
    eTargEEG(:,:,i) = eeg1(:,config.bufferMS+1:(size(eeg1,2)-config.bufferMS));
end 
%sampleRate = GetRateAndFormat(events(1));

 avgTargEEG = squeeze(nanmean(eTargEEG,1));
 avgBaseEEG = squeeze(nanmean(baseEEG,1));
 neTargEEG = NaN(size(avgTargEEG,1),size(avgTargEEG,2));
 for ii=1:size(avgTargEEG,2)
     neTargEEG(:,ii) = (avgTargEEG(:,ii)-nanmean(avgBaseEEG(:,ii)))./std(avgBaseEEG(:,ii));
 end
% 
 avgBackEEG = squeeze(nanmean(eBackEEG,1));
 avgBaseEEG = squeeze(nanmean(baseEEG,1));
 neBackEEG = NaN(size(avgBackEEG,1),size(avgBackEEG,2));
% 
 for ii=1:size(avgBackEEG,2)
     neBackEEG(:,ii) = (avgBackEEG(:,ii)-nanmean(avgBaseEEG(:,ii)))./std(avgBaseEEG(:,ii));
 end
  Early = neTargEEG - neBackEEG;

% 
 b=(0:size(baseEEG,2)-1);
 b=b+ -500;
 g=(0:size(eTargEEG,2)-1);
 g=g+ -500;
 u=(0:size(eBackEEG,2)-1);
 u=u+ -500;
% % 
% %  
 for i = 1:size(eTargEEG, 3)
     for m = 1:size(eTargEEG, 2)
         [H, P, CI, STATS] = ttest2(eTargEEG(:,m,i), eBackEEG(:,m,i));
         t = STATS.tstat ;  
         s = STATS.sd;
         t_mat(m,:,i) = [H;P;CI; t; s]';
 
     end 
 end 

[retain1, retain2] = vplotter(events, 'goLATE');
%if an explant:
%      lateTarg = passEvents(retain1);
%      lateBack = passEvents(retain2);
     %if an implant:
     lateTarg = events(retain1);
     lateBack = events(retain2);     
% 
    lTargEEG = zeros(sum(retain1),config.durationMS,length(elecs)); %%%
    lBackEEG = zeros(sum(retain2),config.durationMS,length(elecs));
    [ev_out]=get_baseline_random(events,config.baseSecBetEv,config.baseJitt); %%%
    %%%change rec_type to lfp but after remember to change it to eeg

for i = 1:length(elecs)
    eeg4 = an_getlfp_ms_wrapper (elecs(i),lateTarg,config.durationMS, config.offsetMS, config.bufferMS, (50.5),'low', 14);
    eeg5 = an_getlfp_ms_wrapper (elecs(i),lateBack,config.durationMS, config.offsetMS, config.bufferMS, (50.5),'low', 14 );
    lTargEEG(:,:,i) = eeg4(:,config.bufferMS+1:(size(eeg4,2)-config.bufferMS));
    lBackEEG(:,:,i) = eeg5(:,config.bufferMS+1:(size(eeg4,2)-config.bufferMS));
end 
%sampleRate = GetRateAndFormat(events(1));

 avgLTargEEG = squeeze(nanmean(lTargEEG,1));
 avgBaseEEG = squeeze(nanmean(baseEEG,1));
 nlTargEEG = NaN(size(avgLTargEEG,1),size(avgLTargEEG,2));
 for ii=1:size(avgLTargEEG,2)
     nlTargEEG(:,ii) = (avgLTargEEG(:,ii)-nanmean(avgBaseEEG(:,ii)))./std(avgBaseEEG(:,ii));
 end
% 
 avgLBackEEG = squeeze(nanmean(eBackEEG,1));
 avgBaseEEG = squeeze(nanmean(baseEEG,1));
 nLBackEEG = NaN(size(avgLBackEEG,1),size(avgLBackEEG,2));
% 
 for ii=1:size(avgLBackEEG,2)
     nLBackEEG(:,ii) = (avgLBackEEG(:,ii)-nanmean(avgBaseEEG(:,ii)))./std(avgBaseEEG(:,ii));
 end
 
 Late = nlTargEEG - nLBackEEG;
% 
 b=(0:size(baseEEG,2)-1);
 b=b+ -500;
 g=(0:size(lTargEEG,2)-1);
 g=g+ -500;
 u=(0:size(lBackEEG,2)-1);
 u=u+ -500;
% % 
% %  
%  for i = 1:size(lTargEEG, 3)
%      for m = 1:size(lTargEEG, 2)
%          [H, P, CI, STATS] = ttest2(lTargEEG(:,m,i), lBackEEG(:,m,i));
%          t2 = STATS.tstat ;  
%          s2 = STATS.sd;
%          t2_mat(m,:,i) = [H;P;CI; t2; s2]';
%  
%      end 
%  end 
%  
saveFigs = 1;
%1 if you want to save, 0 if you want to view
 figDir = fullfile(dirs.scratch,'figs/oddball/neuralynx/',subj,response,'earlyLATE');
 if ~exist(figDir,'dir'),mkdir(figDir);end

for i = 1: size(avgTargEEG,2)
    name= Enames(i,:);
 
    figure ('Name',name); hold on
    plot(g, Early(:,i)); hold on
    plot(u, Late(:,i), 'color', 'm');hold on
%                for r = 1:size(t_mat,1)
%                   q = find(t_mat(:,1,i));
%                   plot (g(q), neTargEEG(q,i), 'y*'); 
%                end
    xlim([-600 700]);
    xlim([-600 700]);
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    line([0,0],ylim, 'color', 'r')
    legend ('Early Trials', 'Late Trials')
      if saveFigs
                  cd(figDir);
                  print(gcf,[ subj '-' name '_goEarlyvLate' ],'-dpng');close

      end
end
