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
feedbackDur = 10;
feedFrames = round(feedbackDur/ifi);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Sync us and get a time stamp
vbl = Screen('Flip', window);
waitframes = 1;

% Maximum priority level
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);

cd /Users/tnl-macbook/Documents/CCDT/data/HUP136/Session_0;
load('sessRTs.mat')
reactTimeMS = sessRTs;
records.avgRT = 350;
records.bestRT = 100;
records.fastCount = 0;
records.subPhysCount = 0;
textSize = 40;

feedbackRTs = reactTimeMS(reactTimeMS(:,1) > 0);
avgRecord = 0; rtRecord = 0;
if round(mean(feedbackRTs)) < records.avgRT
    records.avgRT = round(mean(feedbackRTs));
    avgRecord = 1;
end
if round(min(feedbackRTs)) < records.bestRT
    records.bestRT = round(min(feedbackRTs));
    rtRecord =1;
end
records.fastCount = records.fastCount + sum(feedbackRTs<=250);
records.subPhysCount = records.subPhysCount + sum(feedbackRTs<=150);
save('records.mat','records');

textColors = zeros(4,4);
textColors(:,1) = [round(mean(feedbackRTs)); records.avgRT; round(min(feedbackRTs)); records.bestRT];
for j=1:size(textColors,1)
    if textColors(j,1)<=250
        textColors(j,[2 3 4]) = [0 255 0];
    elseif round(mean(feedbackRTs))<=500
        textColors(j,[2 3 4]) = [255 255 0];
    else
        textColors(j,[2 3 4]) = [255 0 0];
    end
end

for fFrame = 1:feedFrames
    Screen('DrawLine', window, white, screenXpixels*0.5, 0, screenXpixels*0.5, screenYpixels, 2);
    Screen('TextSize',window,80);
    DrawFormattedText(window,'This Set',screenXpixels*0.1,screenYpixels*0.10,white);
    DrawFormattedText(window,'Your Records',screenXpixels*0.6,screenYpixels*0.10,white);
    
    %Average RT this set    %Highest average to date
    Screen('TextSize',window,50);
    DrawFormattedText(window,'Average RT:',screenXpixels*0.1,screenYpixels*0.40,white);
    Screen('TextSize',window,60);
    DrawFormattedText(window,[num2str(round(mean(feedbackRTs))) 'ms'],screenXpixels*0.3,screenYpixels*0.50,textColors(1,[2 3 4]));
    Screen('TextSize',window,50);
    DrawFormattedText(window,'Best Average RT:',screenXpixels*0.6,screenYpixels*0.40,white);
    Screen('TextSize',window,60);
    DrawFormattedText(window,[num2str(records.avgRT) 'ms'],screenXpixels*0.8,screenYpixels*0.50,textColors(2,[2 3 4]));
    if avgRecord
        Screen('TextSize',window,40);
        DrawFormattedText(window,'New record!',screenXpixels*0.6,screenYpixels*0.50,[0 255 0]);
    end
    %Fastest RT this set    %Fastest RT to date
    Screen('TextSize',window,50);
    DrawFormattedText(window,'Fastest RT:',screenXpixels*0.1,screenYpixels*0.70,white);
    Screen('TextSize',window,60);
    DrawFormattedText(window,[num2str(round(min(feedbackRTs))) 'ms'],screenXpixels*0.3,screenYpixels*0.80,textColors(3,[2 3 4]));
    Screen('TextSize',window,50);
    DrawFormattedText(window,'Fastest RT:',screenXpixels*0.6,screenYpixels*0.70,white);
    Screen('TextSize',window,60);
    DrawFormattedText(window,[num2str(records.bestRT) 'ms'],screenXpixels*0.8,screenYpixels*0.80,textColors(4,[2 3 4]));
    if rtRecord
        Screen('TextSize',window,40);
        DrawFormattedText(window,'New record!',screenXpixels*0.6,screenYpixels*0.80,[0 255 0]);
    end
    %Number of fast         %Highest # of fast to date
    %DrawFormattedText(window,fastSetText,screenXpixels*0.1,screenYpixels*0.60,white);
    %DrawFormattedText(window,recordFastText,screenXpixels*0.6,screenYpixels*0.60,white);
    %Number of superhuman   %Highest # of superhuman to date
    %DrawFormattedText(window,subPhysSetText,screenXpixels*0.1,screenYpixels*0.80,white);
    %DrawFormattedText(window,recordSubPhysText,screenXpixels*0.6,screenYpixels*0.80,white);
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
end
%% Keyboard response testing (test usb connected device)
% sca;
% [gpInd, gpName] = GetGamepadIndices;
% %gpName = gpName{1};
% keepGoing = true;
% while keepGoing 
%     buttonPress = Gamepad('GetButton',gpInd,3);
%     [~,~,keyCode] = KbCheck;
%     if buttonPress
%         disp('YAY, that one!')
%         buttonPress = false;
%     end
%     if keyCode(KbName('ESCAPE')), keepGoing = false; end
% end

%% Score Bar testing
% % Score setup
% scoreBarColor = [0 0 1];
% scoreBarWidth = 50;
% %scoreBarBase = [screenXpixels - scoreBarWidth, screenYpixels, screenXpixels, screenYpixels];
% scoreBarBase = [0.05*screenXpixels, 0.75*screenYpixels, 0.05*screenXpixels, 0.75*screenYpixels+scoreBarWidth];
% currentScore = 100;
% trialScore = 1500;
% scoreRect = [0, 0, currentScore, 0];
% scoreBarRect = scoreBarBase + scoreRect;
% disp('score Bar Start Pix')
% disp(0.05*screenXpixels)
% disp('screen x Pix')
% disp(screenXpixels)
% % disp(scoreBarRect)
% if trialScore < 0
%     scoreChangeColor = [1 0 0];
%     scoreChangeRect = [scoreBarRect(3)+trialScore, scoreBarRect(2), scoreBarRect(3), scoreBarRect(4)];
%     newScoreBarRect = [scoreBarRect(1), scoreBarRect(2), scoreChangeRect(1), scoreBarRect(4)];
% else
%     scoreChangeColor = [0 1 0];
%     scoreChangeRect = [scoreBarRect(3), scoreBarRect(2), scoreBarRect(3)+trialScore, scoreBarRect(4)];
%     newScoreBarRect = [scoreBarRect(1), scoreBarRect(2), scoreChangeRect(3), scoreBarRect(4)];
% end
% % disp(scoreChangeRect)
% % disp(newScoreBarRect)
% % Loop the animation until a key is pressed
% for fFrame = 1:feedFrames
%     % Trial report screen with reaction time and initial score bar
%     Screen('TextSize',window,80);
%     DrawFormattedText(window,'250 ms','center','center',scoreBarColor);
%     
%     %Grow/shrink score bar in first half of report screen duration
%     if 2*fFrame < feedFrames
%         if currentScore ~= 0
%             Screen('FillRect', window, scoreBarColor, scoreBarRect);
%         end
%         growRect = scoreChangeRect;
%         if trialScore < 0
%            growRect(1) = scoreBarRect(3)-2*fFrame/feedFrames*(scoreBarRect(3)-scoreChangeRect(1));
%         else
%            growRect(3) = scoreBarRect(3)+2*fFrame/feedFrames*(scoreChangeRect(3)-scoreBarRect(3));
%         end
%         Screen('FillRect', window, scoreChangeColor, growRect);
%     else
%         Screen('FillRect', window, scoreBarColor, newScoreBarRect);
%     end
%     
%     % Flip to the screen
%     vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
% 
% end
%%
% Clear the screen
sca;