function single_trial_power
% Calculate averaged power for individual trials
subj = 'HUP143';
task = 'CCDT';
sessList   = {'Session_0'};
alignment = 'cue'; %cue,CC,response
comparison = 'delay'; %delay, RTcombined, timeCombined(early/late trials), RTlong, RTshort, timeLong, timeShort
delayPeriod = 'both'; %both,long,short
fRange = [3 12]; %freqency band to collapse
tRange = [-500 0];
freqName = sprintf('%d-%d Hz',fRange); %for plot titles
saveFigs = 1;
configNum = 1;
dirs = le_dirs(task);
powConfig = le_config_calcPow(configNum,task);
saveFigsDir = fullfile(dirs.scratch,'figs',subj,'powPlots',num2str(configNum),'singleTrial');

% open jacksheet
jacFile = fullfile(dirs.data,'eeg',subj,'docs','jacksheet.txt');
fid = fopen(jacFile,'r');
JAC=textscan(fid,'%d%s','delimiter','\t');
JAC{:,2} = strtrim(JAC{:,2});
fclose(fid);
channels = JAC{1};
excludeChan = {'EKG1','EKG2','DC1'};
[lia,locb] = ismember(excludeChan,JAC{2});
if sum(lia)>0
    channels(locb(locb~=0)) = [];
end

evFile = fullfile(dirs.events,strcat(subj,'_events'));
load(evFile,'events');
[shortD,longD] = le_events2retain(events,task,alignment,comparison);
switch delayPeriod
    case 'both'
        fixation_events = events(longD|shortD);
    case 'long'
        fixation_events = events(longD);
    case 'short'
        fixation_events = events(shortD);
end
fix_ev_rts = cell2mat({fixation_events.rt});
%add 70-100Hz 500ms precue power and same for 3-12Hz (num_ele x t) array
cat_avg_pow = [];
for chan = 1:length(channels)
    channel = num2str(channels(chan));
    powDir = fullfile(dirs.scratch,'POWER',num2str(configNum),subj,strcat(task,'_events'));
    pow_cat = [];
    for ii=1:length(sessList) %concatenate power for this channel across sessions
        sessPowFile = fullfile(powDir,sessList{ii},channel);
        load(sessPowFile);
        pow_cat = cat(3,pow_cat,pow);
    end
    switch delayPeriod
        case 'both'
            fix_pow = pow_cat(:,:,shortD|longD);
        case 'long'
            fix_pow = pow_cat(:,:,longD);
        case 'short'
            fix_pow = pow_cat(:,:,shortD);
    end
    freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
    [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
    [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
    fInd = fInd_start:fInd_end;
    tBins = nanmean(powConfig.timeBins,2)';
    tf_pow = fix_pow(fInd,tBins>tRange(1)&tBins<tRange(2),:);
    avg_precue_pow = squeeze(nanmean(nanmean(tf_pow,2),1));
    cat_avg_pow = [cat_avg_pow,avg_precue_pow];
end
rt_precue_pow = [fix_ev_rts',nanmean(cat_avg_pow,2)];
sortRt = sort(fix_ev_rts(fix_ev_rts>0));
thresh = [sortRt(round(length(sortRt)/3)),sortRt(round(length(sortRt)/3)*2)];%fastest/slowest 3rd
fastRTs = fix_ev_rts<thresh(1);
slowRTs = fix_ev_rts>thresh(2);
fast_pow = nanmean(cat_avg_pow(fastRTs,:),1);
slow_pow = nanmean(cat_avg_pow(slowRTs,:),1);

figure;
switch delayPeriod
    case 'both'
        title([subj,' All Delay Global Trial Pre-cue Power ',freqName])
    case 'long'
        title([subj,' Long Delay Global Trial Pre-cue Power ',freqName])
    case 'short'
        title([subj,' Short Delay Global Trial Pre-cue Power ',freqName])
end
rt_precue_pow = rt_precue_pow(rt_precue_pow(:,1)>0,:);% filter negative RTs
scatter(rt_precue_pow(:,2),rt_precue_pow(:,1))
xlabel('Power')
ylabel('Reaction Time (ms)')

figure;
switch delayPeriod
    case 'both'
        title([subj,' All Delay Average Pre-cue Power ',freqName])
    case 'long'
        title([subj,' Long Delay Average Pre-cue Power ',freqName])
    case 'short'
        title([subj,' Short Delay Average Pre-cue Power ',freqName])
end
plot(1:size(cat_avg_pow,2),fast_pow,'g')
hold on;
plot(1:size(cat_avg_pow,2),slow_pow,'r')
xlabel('Channel')
ylabel('Power')
legend('Fast','Slow')

keyboard;
