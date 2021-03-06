function CCDTreplay
% This function replays a CCDT session with added gaze position. Reads and
% calibrates eyeData with CCDTcalibrate.m. Press 'esc' to exit at any time.
%% Select subject and session
subj = 'HUP171';
session = 'Session_0';
dataDir = fullfile('/Volumes/HumanStudies/HumanStudies/CCDT/eeg',subj);
taskImgDir = '/Volumes/HumanStudies/MATLAB/LucasLabHumanStudies/MATLAB/CCDT';
suppressGazeDisp = 1; % only a few subjects have real gaze data
liveDataFeed = 1;
cgaze = [255 0 0]; % gaze dot color
gazeSize = 10; % gaze dot size
load(fullfile(dataDir,'behavioral',session,'sessConfig.mat'));
load(fullfile(dataDir,'behavioral',session,'sessRTs.mat'));
load(fullfile(dataDir,'behavioral',session,'targetPix.mat'));

%% color info
sca;
PsychDefaultSetup(2);
screenNumber = max(Screen('Screens')); % use 2nd monitor if present
white = WhiteIndex(screenNumber); % color of fixation/response crosses
black = BlackIndex(screenNumber); % color of screen
cfix = ntsc2rgb([0.8 0 0]); %sessConfig.cfix % color of fixation target
cres = ntsc2rgb([0.8 .3 -.3]); %sessConfig.cres % color of response target
% colors for feedback text
cfast = sessConfig.cfast; % green
cmid = sessConfig.cmid; % yellow
cslow = sessConfig.cslow; % red

%% position info
partialScreen = [0,0,400,400]; % Window size for debugging
Screen('Preference','SkipSyncTests', 1); % don't worry about accurate timing
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black, partialScreen);
[xCenter, yCenter] = RectCenter(windowRect);
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
screenScale = mean([screenXpixels/1440,screenYpixels/900]);
fixDim = round(sessConfig.fixDim*screenScale); % fixation target diameter (pixels)
baseRect = [0 0 fixDim fixDim];
textSize = round(sessConfig.textSize*screenScale);

if screenScale ~= 1,targetPix = round(targetPix*screenScale); end % adjust targetPix for new resolution
fixCoord = NaN*ones(4,length(targetPix));
for ii = 1:length(targetPix)
    fixCoord(:,ii) = CenterRectOnPointd(baseRect, targetPix(ii,1), targetPix(ii,2));
end
subPhysImage = imread([taskImgDir,'/subPhysFeedback.jpg']);
scaledSubPhysRect = [0 0 0.5*size(subPhysImage,2) 0.5*size(subPhysImage,1)];
subPhysRect = CenterRectOnPointd(scaledSubPhysRect,screenXpixels/2,0.5*size(subPhysImage,1));
subPhysTexture = Screen('MakeTexture', window, subPhysImage);

%% Calibrate and condition gaze data
[gazeP,eeg,fs] = CCDTcalibrate(dataDir,session,liveDataFeed);
sessMS = size(eeg,2)/fs*1000;
totFrames = sum(~ismember(gazeP(:,3),[98 99])); %some lines don't correspond to frame changes
rr = sessMS/totFrames; %screen refresh rate = ms per gaze sample (avg for session)
rrS = fs/1000*rr; %samples per screen refresh
[b,a] = butter(2,10/(60/2));
gazeP(:,1:2) = filtfilt(b,a,gazeP(:,1:2));

%% Data analysis setup
%trial window settings
durationMS=3000;
offsetMS=-1000;
bufferMS=1000;
preCueSamples = fs/1000*abs(offsetMS);
postCueSamples = fs/1000*durationMS+offsetMS;
bufferSamples = fs/1000*bufferMS;

%% timing info
ifi = Screen('GetFlipInterval', window); % inter-frame interval (s)
waitframes = 1; % numer of frames to wait before re-drawing
rspDurC = sessConfig.rspDurC; % response interval duration - CCDT (sec)
subPhysFeedB = sessConfig.subPhysFeedB;
fastFeedB = sessConfig.fastFeedB;
slowFeedB = sessConfig.slowFeedB;

%% score setup
scoreBarColor = sessConfig.scoreBarColor;
scoreBarWidth = sessConfig.scoreBarWidth;
scoreBarBase = [0.05*screenXpixels, 0.75*screenYpixels, 0.05*screenXpixels, 0.75*screenYpixels+scoreBarWidth];
currentScore = 0;
barAvailable = false; % don't call to display if score = 0. will throw error
scoreBase = sessConfig.scoreBase;
slowScore = sessConfig.slowScore;

%% Experimental loop
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);
HideCursor();
vbl = Screen('Flip', window);
frame = 1;
bf=0;% bad frame counter (doesn't correspond to actual frame change)
target = 1; % target = index of fixFrame
for i=1: length(sessRTs)
    %% update data analysis with previous trial
    if i>1,CCDT_train(eeg(:,cueInd-preCueSamples:cueInd+postCueSamples),fs,bufferSamples,cueInd,respInd,sessRTs(i-1));end 
    
    if gazeP(frame,3) == 99, frame = frame+1;bf=bf+1; end % Start of loop, does not correspond to a frame
    if gazeP(frame,3) == 98, frame = frame+1;bf=bf+1; end % Start of trial, does not correspond to a frame
    scoreRect = [0, 0, currentScore, 0];
    scoreBarRect = scoreBarBase + scoreRect;
    %% inter trial interval
    while gazeP(frame,3) == 0
        if ~suppressGazeDisp
            if gazeP(frame,1)>0,Screen('DrawDots', window, gazeP(frame,1:2),gazeSize,cgaze,[],1);end%display gaze
        end
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        frame = frame+1;
        [~,~,keyCode] = KbCheck;
        if keyCode(KbName('ESCAPE')),sca; return; end
    end
    %% fixation interval
    loopTarg = gazeP(frame,3); %target # within loop (i.e. 1-Ntarg)
    cueInd = round(rrS*(frame-bf));
    while gazeP(frame,3) == loopTarg
        Screen('FillRect', window, cfix, fixCoord(:,target));
        if ~suppressGazeDisp
            if gazeP(frame,1)>0,Screen('DrawDots', window, gazeP(frame,1:2),gazeSize,cgaze,[],1);end%display gaze
        end
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        frame = frame+1;
        [~,~,keyCode] = KbCheck;
        if keyCode(KbName('ESCAPE')), sca; return; end
    end
    respInd = round(rrS*(frame-bf));
    if gazeP(frame,3)==98,continue;end % early press, repeat target
    
    %% response interval
    while gazeP(frame,3) > loopTarg
        Screen('FillRect', window, cres, fixCoord(:,target));
        if ~suppressGazeDisp
            if gazeP(frame,1)>0,Screen('DrawDots', window, gazeP(frame,1:2),gazeSize,cgaze,[],1);end%display gaze
        end
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        frame = frame+1;
        [~,~,keyCode] = KbCheck;
        if keyCode(KbName('ESCAPE')),sca; return; end
    end
    if sessRTs(i,1)>0 && sessRTs(i,1) < rspDurC % Skip feedback on early press/non-response
        % Set feedback info
        trialRT = sessRTs(i,1);
        trialScore = round((10*(1-trialRT/500)*scoreBase)*screenScale);
        if trialRT <= subPhysFeedB
            feedbackcolor = cfast;
        elseif trialRT <= fastFeedB
            feedbackcolor = cfast;
        elseif trialRT<= slowFeedB
            feedbackcolor = cmid;
        else
            feedbackcolor = cslow;
            trialScore = round(slowScore*screenScale);
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
        %% feedback interval
        growStart = frame+1;
        growStop = frame + (find(gazeP(frame:end,3),1))/2;
        if isempty(growStop),growStop = length(gazeP);end 
        while gazeP(frame,3)==0
            Screen('TextSize',window,textSize);
            if trialRT <= subPhysFeedB,Screen('DrawTexture', window, subPhysTexture,[],subPhysRect,0);end
            DrawFormattedText(window,[num2str(round(trialRT)) 'ms'],'center','center',feedbackcolor);
            if frame < growStop && frame > growStart
                if barAvailable
                    Screen('FillRect', window, scoreBarColor, scoreBarRect);
                end
                growRect = scoreChangeRect;
                if trialScore < 0
                    growRect(1) = scoreBarRect(3)-(frame-growStart)/(growStop-growStart)*(scoreBarRect(3)-scoreChangeRect(1));
                    if barAvailable, Screen('FillRect', window, feedbackcolor, growRect);end
                elseif trialScore > 0
                    growRect(3) = scoreBarRect(3)+(frame-growStart)/(growStop-growStart)*(scoreChangeRect(3)-scoreBarRect(3));
                    if barAvailable, Screen('FillRect', window, feedbackcolor, growRect);end
                end
            else
                if barAvailable, Screen('FillRect', window, scoreBarColor, newScoreBarRect);end
            end
            if ~suppressGazeDisp
                if gazeP(frame,1)>0,Screen('DrawDots', window, gazeP(frame,1:2),gazeSize,cgaze,[],1);end%display gaze
            end
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            frame = frame+1;
            if frame > length(gazeP), return; end
            [~,~,keyCode] = KbCheck;
            if keyCode(KbName('ESCAPE')), sca; return; end
        end
        target = target+1;
    end
end
Priority(0);
sca;
end