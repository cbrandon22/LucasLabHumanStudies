function CCDT  
% Color-change detection task
% Set subject/session info for correct directories and select delay paradigm.
% Before beginning, task will enter loop while sending sync pulses. Task
% begins after pressing any key. Press 'esc' to end task.
% Sync pulses/eye tracker dependent on labjack.m
% Task will save:
%       sessRTs(reaction time,response time, fixation duration)
%       syncTimes(time of each TTL pulse)
%       eyeData(xVolts,yVolts,diameterVolts,targetNumber)
%       targetPix(xPixels,yPixels) index = targetNumber
% Saves every 9 successful trials

%%    Subject Info
% Set subject and session, if continuing session, load variables
subj = 'test3';
session = 'Session_0';
delayParadigm = 'discrete'; %Set to 'fixed','discrete',or 'continuous'
Ncal = 6; % Ncal = repetitions per calibration target

sessdir = fullfile('~/Documents/CCDT/data',subj,session);
if exist(sessdir,'dir') ~=7 %check if first entry in session
    sessRTs = []; %initialize session reation time array
    syncTimes = []; %initialize sync times array
    eyeData = []; %initialize eye tracking voltage array
    targetPix = [];
elseif exist(sessdir,'dir') ==7 %if session has already been started, load variables
    cd (sessdir);
    load('sessRTs.mat');
    load('syncTimes.mat');
    load('targetPix.mat');
    load('eyeData.mat');
end

%% Color info
sca;
PsychDefaultSetup(2);
screenNumber = max(Screen('Screens')); % use 2nd monitor if present
white = WhiteIndex(screenNumber); % color of fixation/response crosses
black = BlackIndex(screenNumber); % color of screen
cfix = ntsc2rgb([0.8 0 0]); % color of fixation target
cres = ntsc2rgb([0.8 .3 -.3]); % color of response target
% colors for feedback text
cfast = [0 255 0];% green
cmid = [255 255 0];% yellow
cslow = [255 0 0]; % red

%% LabJack setup
% Creates labjack object using labjack.m. Set silentMode to true in object
% construction ('silentMode',true) to run without communicating with LabJack or saving TTL/eye
% tracker data. Automatically enters silentMode if no LabJack connected.
daq = labJack('verbose',true);

% LabJack FIO port settings
daqinfo.pulsePort = 0; %TTL pulse port output
daqinfo.xPort = 4; %x coordinate voltage input
daqinfo.yPort = 5; %y coordinate voltage input
daqinfo.diamPort = 6; %pupil diameter voltage input
daqinfo.pulsewidth = 5; %milliseconds
daqinfo.recBuff = zeros(16,1); % Read bytes (length = 10+ 2*Num of channels to read AIN)

% ConfigIO byte. See Labjack U3 documentation (5.2.3)
% Sets all channels other than FIO0 and FIO1 to analog input, Use ioBuffer[11] to change channel
% config (8 bit command, 0=digital,1=analog ex: 252 = 11111100)
ioBuffer = [0 248 3 11 0 0 13 0 64 0 252 0]';
ioBuffer = daq.checksum(ioBuffer,'extended');
bytesWritten = daq.rawWrite(ioBuffer);
if bytesWritten == 0, warning('ConfigIO error: No bytes written.');end
if bytesWritten < 12, warning('ConfigIO error: Did not write all of the buffer.');end
recBuffer = zeros(12,1);
[bytesRead,retBytes] = daq.rawRead(recBuffer,12);
if bytesRead == 0, warning('ConfigIO error: Read failed.');end
if bytesRead < 12, warning('ConfigIO error: Did not read all of the buffer.');end

% Bytes to send to read AIN from eye tracker channels
sendBuff = [0 248 5 0 0 0 0 1 daqinfo.xPort 31 1 daqinfo.yPort 31 1 daqinfo.diamPort 31]';
sendBuff = daq.checksum(sendBuff,'extended');

%% Position info
fixDim = 16;        % fixation target diameter (pixels)
picDim = 267;       % inter-target length (pixels)
partialScreen = [0,0,400,400]; % Window size for debugging
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);% add partialScreen to debug
[xCenter, yCenter] = RectCenter(windowRect);
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
[Px,Py] = meshgrid(-1:1,-1:1); Ntarg = 9; % 9-point calibration
Px = Px(:)*picDim+xCenter; Py = Py(:)*picDim+yCenter;
baseRect = [0 0 fixDim fixDim];
fixCoord = NaN*ones(4,Ntarg);
for ii = 1:Ntarg
    fixCoord(:,ii) = CenterRectOnPointd(baseRect, Px(ii), Py(ii));
end
picRect = [0 0 2*picDim 2*picDim];
picCoord = CenterRectOnPointd(picRect, screenXpixels/2, screenYpixels/2);
cd ~/Documents/MATLAB/CCDT;
subPhysImage = imread('subPhysFeedback.jpg');
scaledSubPhysRect = [0 0 0.5*size(subPhysImage,2) 0.5*size(subPhysImage,1)];
subPhysRect = CenterRectOnPointd(scaledSubPhysRect,screenXpixels/2,0.5*size(subPhysImage,1));

%% Timing info
switch lower(delayParadigm)
    case 'fixed'
        fixRng = [1.0];
        iscontinuous = false;
    case 'discrete'
        fixRng = [0.5 1.0 1.5];
        iscontinuous = false;
    case 'continuous'
        fixRng = [0.5 1.5];
        isContinuous = true;
end
rspDurC = 1;        % response interval duration - CCDT (sec)
feedbackDur = 2;    % feedback display duration (sec)
stimTestInt = 1;    % Pulse interval during "Press any Key" screen
iti = [2.5 3.5];          % intertrial interval duration (sec)
ifi = Screen('GetFlipInterval', window); % inter-frame interval (s)
rspFramesC = round(rspDurC/ifi);
feedFrames = round(feedbackDur/ifi);
stimFrames = round(stimTestInt/ifi);
waitframes = 1; % number of frames to wait before re-drawing
subPhysFeedB = 150;
fastFeedB = 250;
slowFeedB = 500;

%% Score setup
scoreBarColor = [0 0 1];
scoreBarWidth = 50;
scoreBarBase = [0.05*screenXpixels, 0.75*screenYpixels, 0.05*screenXpixels, 0.75*screenYpixels+scoreBarWidth];
currentScore = 0;
barAvailable = false; % don't call to display if score = 0. will throw error
subPhysScore = 32;
fastScore = 16;
midScore = 8;
slowScore = -10;

Ntpr = Ntarg; % Ntpr = targets per repetition
totalRT = 0;

Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);
vbl = Screen('Flip', window);
subPhysTexture = Screen('MakeTexture', window, subPhysImage);

%% Loop while sending pulses every 1s until keypress
keypress = false;
while ~keypress 
    for frame=1:stimFrames
        Screen('TextSize',window,80);
        DrawFormattedText(window,'Press any key to begin','center','center',white);
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        [keypress,~,keyCode] = KbCheck;
        if keyCode(KbName('ESCAPE')), sca; return; end
        if keypress, break;end
    end
    sendPulse;
end
HideCursor();

%% experimental loop
for irep = 1:Ncal
    targ = randperm(Ntarg);
    itarg = 1;
    reactTimeMS = []; %initialize reation time array
    targetPix = [targetPix;[Px(targ),Py(targ)]];
    readInput(99);% Set target to 99 to indicate loop start

    while itarg<=min(Ntpr,Ntarg)
        respbit = 0; % reset response bit
        readInput(98);% Set target to 98 to indicate trial start
        
        % Create ScoreBarRect with current score
        scoreRect = [0, 0, currentScore, 0];
        scoreBarRect = scoreBarBase + scoreRect;
        
        % intertrial interval
        itiDur = iti(1) + diff(iti)*rand;
        itiFrames = round(itiDur/ifi);
        for frame = 1:itiFrames
            if frame ==1, sendPulse; end
            readInput(0);
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            [~,~,keyCode] = KbCheck;
            if keyCode(KbName('ESCAPE')), respbit = 1; sca; return; end
        end
        
        % fixation interval
        if iscontinuous
            fixDur = fixRng(1) + diff(fixRng)*rand;
        else
            fixDur = datasample(fixRng,1);
        end
        fixFrames = round(fixDur/ifi);
        for frame = 1:fixFrames
            Screen('FillRect', window, cfix, fixCoord(:,targ(itarg)));
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            readInput(itarg);
            if frame == 1, trialStart = vbl; end% get stim onset in ms
            if frame<4 && rem(frame,2) ~=0, sendPulse; end %1,3
            % Read response bit
            [~,rspTime,keyCode] = KbCheck;
            if keyCode(KbName('DownArrow'))
                respbit = 1;
                stimMS = (trialStart + fixDur)*1000;
                trialRT = rspTime*1000-stimMS;
                trialPressTime = rspTime*1000;
                trialDelayTime = fixDur*1000;
                break; % Restart trial if early press
            end
            if keyCode(KbName('ESCAPE')), respbit = 1; sca; return; end
        end
        if respbit
            reactTimeMS = [reactTimeMS;[trialRT,trialPressTime,trialDelayTime]];
            continue;
        end
        
        % response interval
        for frame = 1:rspFramesC
            readInput(itarg+0.1);
            if frame<6 && rem(frame,2) ~=0, sendPulse; end %1,3,5
            Screen('FillRect', window, cres, fixCoord(:,targ(itarg)));
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            % Read response bit
            [~,rspTime,keyCode] = KbCheck;
            if keyCode(KbName('DownArrow'))
                respbit = 1;
                stimMS = (trialStart + fixDur)*1000;
                trialRT = rspTime*1000-stimMS;
                trialPressTime = rspTime*1000;
                trialDelayTime = fixDur*1000;
                totalRT = totalRT + trialRT;
                % Set feedback info
                if trialRT <= subPhysFeedB
                    feedbackcolor = cfast;
                    trialScore = subPhysScore;
                elseif trialRT >subPhysFeedB && trialRT <= fastFeedB
                    feedbackcolor = cfast;
                    trialScore = fastScore;
                elseif trialRT >fastFeedB && trialRT<= slowFeedB
                    feedbackcolor = cmid;
                    trialScore = midScore;
                elseif trialRT > slowFeedB
                    feedbackcolor = cslow;
                    trialScore = slowScore;
                    scoreChangeRect = [scoreBarRect(3)+trialScore, scoreBarRect(2), scoreBarRect(3), scoreBarRect(4)];
                    newScoreBarRect = [scoreBarRect(1), scoreBarRect(2), scoreChangeRect(1), scoreBarRect(4)];
                end
                if trialScore > 0
                    barAvailable = true;
                    scoreChangeRect = [scoreBarRect(3), scoreBarRect(2), scoreBarRect(3)+trialScore, scoreBarRect(4)];
                    newScoreBarRect = [scoreBarRect(1), scoreBarRect(2), scoreChangeRect(3), scoreBarRect(4)];
                end
                currentScore = currentScore + trialScore; % Recalculate startScore
                if currentScore < 0 % Do not set negative startScore
                    currentScore = 0;
                    barAvailable = false; %reset to unavailable
                end
                for fFrame = 1:feedFrames
                    readInput(0);
                    Screen('TextSize',window,80);
                    if trialRT <= subPhysFeedB,Screen('DrawTexture', window, subPhysTexture,[],subPhysRect,0);end
                    DrawFormattedText(window,[num2str(round(trialRT)) 'ms'],'center','center',feedbackcolor);
                    if fFrame > 5 && 2*fFrame < feedFrames
                        if barAvailable
                            Screen('FillRect', window, scoreBarColor, scoreBarRect);
                        end
                        growRect = scoreChangeRect;
                        if trialScore < 0
                            growRect(1) = scoreBarRect(3)-2*fFrame/feedFrames*(scoreBarRect(3)-scoreChangeRect(1));
                        else
                            growRect(3) = scoreBarRect(3)+2*fFrame/feedFrames*(scoreChangeRect(3)-scoreBarRect(3));
                        end
                        if barAvailable, Screen('FillRect', window, feedbackcolor, growRect);end
                    else
                        if barAvailable, Screen('FillRect', window, scoreBarColor, newScoreBarRect);end
                    end
                    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                end
                break;
            end
            if keyCode(KbName('ESCAPE')), respbit = 1; sca; return; end
        end
        if respbit
            reactTimeMS = [reactTimeMS;[trialRT,trialPressTime,trialDelayTime]];
            itarg = itarg + 1;
        end
    end
    if exist(sessdir,'dir') ~=7 %check if first entry in session
        mkdir(sessdir);
    end
    cd (sessdir);
    sessRTs = [sessRTs;reactTimeMS];
    save('sessRTs.mat','sessRTs');
    save('syncTimes.mat','syncTimes');
    save('targetPix.mat','targetPix');
    save('eyeData.mat','eyeData');
end
%% Summary Feedback
if mean(reactTimeMS(:,1)) <= subPhysFeedB
    feedbackcolor = cfast;
elseif mean(reactTimeMS(:,1)) <= fastFeedB
    feedbackcolor = cfast;
elseif mean(reactTimeMS(:,1))<= slowFeedB
    feedbackcolor = cmid;
else
    feedbackcolor = cslow;
end
reactionTimes = sessRTs(:,1);
posRTs = reactionTimes(reactionTimes>0);
finalText = ['Average Reaction Time: ' num2str(round(mean(posRTs))) 'ms\n\n Total Score:  ' num2str(currentScore) ' points!' '\n\n\n\npress any key to exit'];
Screen('TextSize',window,80);
DrawFormattedText(window,finalText,'center','center',feedbackcolor);
Screen('Flip',window);
KbStrokeWait;
Priority(0);
sca;

%% Local functions
function sendPulse
    % Send TTL pulse
    if daq.silentMode, return; end
    time = GetSecs; %getting timestamp. Cmd-response for labjack is <1ms
    daq.timedTTL(daqinfo.pulsePort, daqinfo.pulsewidth);
    syncTimes = [syncTimes; time];
end

function readInput(target)
    % Writes to then reads from labjack for eyetracker input
    if daq.silentMode, return; end
    bytesWritten = daq.rawWrite(sendBuff);
    if bytesWritten == 0, warning('Feedback error: No bytes written.');end
    if bytesWritten < 16, warning('Feedback error: Did not write all of the buffer.');end
    [bytesRead, retBytes] = daq.rawRead(daqinfo.recBuff,16);
    if bytesRead == 0, warning('Feedback error: Read failed.');end
    if bytesRead < length(daqinfo.recBuff), warning('Feedback error: Did not read all of the buffer.');end
    xVolts = double(retBytes(11))*256 + double(retBytes(10));
    yVolts = double(retBytes(13))*256 + double(retBytes(12));
    diamVolts = double(retBytes(15))*256 + double(retBytes(14));
    eyeData = [eyeData;[xVolts,yVolts,diamVolts,target]];
end
end