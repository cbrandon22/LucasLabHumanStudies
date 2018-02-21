function CCDTbehavioral
% Quick look at behavioral data for CCDT
subj = 'HUP133';
sessNum = [0]; % Pick a specific session. Default combines all sessions
dirs = le_dirs;
dataDir = fullfile(dirs.data,'eeg','CCDT',subj,'behavioral');

if isempty(sessNum)
    sessNum = 0;
    allRTs = [];
    sessInd = [1];
    while exist(fullfile(dataDir,['Session_' num2str(sessNum)]),'dir')==7
        cd(fullfile(dataDir,['Session_' num2str(sessNum)]));
        load('sessRTs.mat');
        allRTs = [allRTs; sessRTs];
        sessInd = [sessInd length(allRTs)+1];
        sessNum = sessNum + 1;
    end
else
    cd(fullfile(dataDir,['Session_' num2str(sessNum)]));
    load('sessRTs.mat');
    allRTs = sessRTs;
end

reactionTimes = allRTs(:,1); % only reaction times
posRTs = reactionTimes(reactionTimes>0); % exclude early presses
correctResponseRTs = reactionTimes(reactionTimes>-150); % exclude reactions to fixation start

%Scatter plot of all RT's with best fit line
trials = 1:length(reactionTimes);
scatter(trials,reactionTimes);
hold all;
pcoef = polyfit(trials',reactionTimes,1);
plot(polyval(pcoef,trials));
title('HUP133 Fixed Interval, All RTs');
xlabel('trial number');
ylabel('reaction time (ms)');

%Scatter plot of all RT's with best fit line
trials = 1:length(correctResponseRTs);
scatter(trials,correctResponseRTs);
hold all;
pcoef = polyfit(trials',correctResponseRTs,1);
plot(polyval(pcoef,trials));
title('HUP133 Fixed Interval, "correct" RTs');
xlabel('trial number');
ylabel('reaction time (ms)');

% Average RT per 9-target set
i=9;
sessionAverages = [];
while i<=length(posRTs)
    thisAvg = mean(posRTs((i-8):i));
    sessionAverages = [sessionAverages;thisAvg];
    i=i+9;
end
plot(sessionAverages)
title('HUP133 Fixed Interval, positive average RTs');
xlabel('9-target set number');
ylabel('average reaction time (ms)');

% Windowed average
windowSize = 5;
stepSize = 3;
i = windowSize;
movingAvg = [];
while i<=length(posRTs)
    winAvg = mean(posRTs(i-(windowSize-1):i));
    movingAvg = [movingAvg;winAvg];
    i = i+stepSize;
end
plot(movingAvg)

% Number of guesses per 9-target set (sub 150ms RT)
sessGuess = [];
goodTrial = 0;
sessionGuess = 0;
for i=1:length(allRTs)
    if allRTs(i,1)<150 && allRTs(i,1) >-150
        sessionGuess = sessionGuess+1;
    end
    if allRTs(i,1)>0
        goodTrial = goodTrial+1;
    end
    if rem(goodTrial,9)==0
        sessGuess = [sessGuess;sessionGuess];
        goodTrial = 0;
        sessionGuess = 0;
    end
end
x = [1:length(sessGuess)];
scatter(x,sessGuess)
title('HUP-LE Fixed Interval, guesses per 9-target set');
xlabel('9-trial set number');
ylabel('number of guesses');

% count fast-RTs per 9-target set
sessFast = [];
goodTrial = 0;
sessionFast = 0;
for i=1:length(allRTs)
    if allRTs(i,1)<250 && allRTs(i,1) >150
        sessionFast = sessionFast+1;
    end
    if allRTs(i,1)>0
        goodTrial = goodTrial+1;
    end
    if rem(goodTrial,9)==0
        sessFast = [sessFast;sessionFast];
        goodTrial = 0;
        sessionFast = 0;
    end
end
x = [1:length(sessFast)];
scatter(x,sessFast)
title('HUP-LE Fixed Interval, Fast-RTs per 9-target set');
xlabel('9-trial set number');
ylabel('number of fast-RTs');

% count mid-RTs per 9-target set
sessMid = [];
goodTrial = 0;
sessionMid = 0;
for i=1:length(sessRTs)
    if sessRTs(i,1)<500 && sessRTs(i,1) >250
        sessionMid = sessionMid+1;
    end
    if sessRTs(i,1)>0
        goodTrial = goodTrial+1;
    end
    if rem(goodTrial,9)==0
        sessMid = [sessMid;sessionMid];
        goodTrial = 0;
        sessionMid = 0;
    end
end
x = [1:length(sessMid)];
scatter(x,sessMid)
title('HUP-LE Fixed Interval, Fast-RTs per 9-target set');
xlabel('9-trial set number');
ylabel('number of fast-RTs');
end