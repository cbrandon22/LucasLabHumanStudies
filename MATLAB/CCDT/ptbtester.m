function ptbtester

% Clear the workspace and the screen
sca;
close all;
clearvars;

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);
feedbackDur = 2;
feedFrames = round(feedbackDur/ifi);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Our square will oscilate with a sine wave function to the left and right
% of the screen. These are the parameters for the sine wave
% See: http://en.wikipedia.org/wiki/Sine_wave
amplitude = screenXpixels * 0.25;
frequency = 0.2;
angFreq = 2 * pi * frequency;
startPhase = 0;
time = 0;

% Sync us and get a time stamp
vbl = Screen('Flip', window);
waitframes = 1;

% Maximum priority level
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);

% Score setup
scoreBarColor = [0 0 1];
scoreBarWidth = 50;
%scoreBarBase = [screenXpixels - scoreBarWidth, screenYpixels, screenXpixels, screenYpixels];
scoreBarBase = [0.05*screenXpixels, 0.75*screenYpixels, 0.05*screenXpixels, 0.75*screenYpixels+scoreBarWidth];
currentScore = 100;
trialScore = -20;
scoreRect = [0, 0, currentScore, 0];
scoreBarRect = scoreBarBase + scoreRect;
disp(scoreBarRect)
if trialScore < 0
    scoreChangeColor = [1 0 0];
    scoreChangeRect = [scoreBarRect(3)+trialScore, scoreBarRect(2), scoreBarRect(3), scoreBarRect(4)];
    newScoreBarRect = [scoreBarRect(1), scoreBarRect(2), scoreChangeRect(1), scoreBarRect(4)];
else
    scoreChangeColor = [0 1 0];
    scoreChangeRect = [scoreBarRect(3), scoreBarRect(2), scoreBarRect(3)+trialScore, scoreBarRect(4)];
    newScoreBarRect = [scoreBarRect(1), scoreBarRect(2), scoreChangeRect(3), scoreBarRect(4)];
end
disp(scoreChangeRect)
disp(newScoreBarRect)
% Loop the animation until a key is pressed
for fFrame = 1:feedFrames
    % Trial report screen with reaction time and initial score bar
    Screen('TextSize',window,80);
    DrawFormattedText(window,'250 ms','center','center',scoreBarColor);
    
    %Grow/shrink score bar in first half of report screen duration
    if 2*fFrame < feedFrames
        if currentScore ~= 0
            Screen('FillRect', window, scoreBarColor, scoreBarRect);
        end
        growRect = scoreChangeRect;
        if trialScore < 0
           growRect(1) = scoreBarRect(3)-2*fFrame/feedFrames*(scoreBarRect(3)-scoreChangeRect(1));
        else
           growRect(3) = scoreBarRect(3)+2*fFrame/feedFrames*(scoreChangeRect(3)-scoreBarRect(3));
        end
        Screen('FillRect', window, scoreChangeColor, growRect);
    else
        Screen('FillRect', window, scoreBarColor, newScoreBarRect);
    end
    
%     % Position of the square on this frame
%     xpos = amplitude * sin(angFreq * time + startPhase);
% 
%     % Add this position to the screen center coordinate. This is the point
%     % we want our square to oscillate around
%     squareXpos = xCenter + xpos;
% 
%     % Center the rectangle on the centre of the screen
%     centeredRect = CenterRectOnPointd(baseRect, squareXpos, yCenter);
% 
%     % Draw the rect to the screen
%     Screen('FillRect', window, rectColor, centeredRect);
 
    % Flip to the screen
    vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);

end

% Clear the screen
sca;