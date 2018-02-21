function CCDTreplay(subj,session)
dataDir = fullfile('~/Documents/CCDT/data',subj,session);
cd(dataDir);
load('sessRTs');
load('targetPix');
gazeP = CCDTcalibrate(subj,session);

% color info
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

% position info
fixDim = 16;        % fixation target diameter (pixels)
picDim = 267;       % inter-target length (pixels)
window = PsychImaging('OpenWindow', screenNumber, black);
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
baseRect = [0 0 fixDim fixDim];
fixCoord = NaN*ones(4,length(targetPix));
for ii = 1:length(targetPix)
    fixCoord(:,ii) = CenterRectOnPointd(baseRect, targetPix(ii,1), targetPix(ii,2));
end
cd ~/Documents/MATLAB/CCDT;
subPhysImage = imread('subPhysFeedback.jpg');
scaledSubPhysRect = [0 0 0.5*size(subPhysImage,2) 0.5*size(subPhysImage,1)];
subPhysRect = CenterRectOnPointd(scaledSubPhysRect,screenXpixels/2,0.5*size(subPhysImage,1));
subPhysTexture = Screen('MakeTexture', window, subPhysImage);

% timing info
ifi = Screen('GetFlipInterval', window); % inter-frame interval (s)
waitframes = 1; % numer of frames to wait before re-drawing
subPhysFeedB = 150;
fastFeedB = 250;
slowFeedB = 500;

% score setup
scoreBarColor = [0 0 1];
scoreBarWidth = 50;
scoreBarBase = [0.05*screenXpixels, 0.75*screenYpixels, 0.05*screenXpixels, 0.75*screenYpixels+scoreBarWidth];
currentScore = 0;
barAvailable = false; % don't call to display if score = 0. will throw error
subPhysScore = 32;
fastScore = 16;
midScore = 8;
slowScore = -10;

% Experimental loop
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);
vbl = Screen('Flip', window);
frame = 1;
for i=1: length(sessRTs)
    if gazeP(frame,3) == 99, frame = frame+1; end % Start of loop, does not correspond to a frame
    if gazeP(frame,3) == 98, frame = frame+1; end % Start of trial, does not correspond to a frame
    %inter trial interval
    while gazeP(frame,3) == 0
        %display gaze
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        frame = frame+1;
    end
    %fixation interval
    target = gazeP(frame,3);
    while gazeP(frame,3) == target
        %display gaze
        Screen('FillRect', window, cfix, fixCoord(:,target));
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        frame = frame+1;
    end
    %response interval
    while gazeP(frame,3) > target
        %display gaze
        Screen('FillRect', window, cres, fixCoord(:,targ(itarg)));
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        frame = frame+1;
    end
    % Set feedback info
    trialRT = sessRTs(i,1);
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
    %feedback interval
    barStart = frame+5;
    barStop = frame + find(gazeP(frame:end,3),1)-2;
    while gazeP(frame,3)==0
        %display gaze
        Screen('TextSize',window,80);
        if trialRT <= subPhysFeedB,Screen('DrawTexture', window, subPhysTexture,[],subPhysRect,0);end
        DrawFormattedText(window,[num2str(round(trialRT)) 'ms'],'center','center',feedbackcolor);
        if frame < barStop && frame > barStart
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
        frame = frame+1;
    end
end
Priority(0);
sca;
end