function CCDT_noFeedback  
% Color-change detection task
% Set subject/session info for correct directories and select delay paradigm.
% Before beginning, task will enter loop while sending sync pulses. Task
% begins after pressing any key. Press 'esc' to end task.
% Sync pulses/eye tracker dependent on labjack.m
% Task will save:
%       sessConfig(variable settings for CCDTreplay)
%       sessRTs(reaction time,response time, fixation duration)
%       syncTimes(time of each TTL pulse)
%       eyeData(xVolts,yVolts,diameterVolts,targetNumber)
%       targetPix(xPixels,yPixels) index = targetNumber
% Saves every 9 successful trials
% Written by Cameron Brandon 10/2016 (cbrandon22@gmail.com)

%% Subject Info
% Set subject and session, if continuing session, load variables
subj = 'testScore2';
session = 'Session_0';
delayParadigm = 'discrete'; %Set to 'fixed','discrete',or 'continuous'
Ncal = 6; % Ncal = repetitions per calibration target

recordsDir = fullfile('~/Documents/CCDT/data',subj); % Where to store subj's records for feedback
if exist(recordsDir,'dir') ~= 7,mkdir(recordsDir);end
cd(recordsDir);
if exist('records.mat','file') ~=2 %initialize or load subj results
    records.avgRT = 1000;
    records.bestRT = 1000;
    records.fastCount =0;
    records.subPhysCount = 0;
    save('records.mat','records');
else
    cd(recordsDir);
    load('records.mat');
end
sessdir = fullfile('~/Documents/CCDT/data',subj,session);
if exist(sessdir,'dir') ~=7 %check if first entry in session
    sessRTs = []; %session reation times, key press times, delay periods
    syncTimes = []; %sync times
    eyeData = []; %eye tracking voltages
    targetPix = []; %target pixel coordinates
elseif exist(sessdir,'dir') ==7 %if session has already been started, load variables
    cd (sessdir);
    load('sessRTs.mat');
    load('syncTimes.mat');
    load('targetPix.mat');
    load('eyeData.mat');
end

%% Color info
%sca;
PsychDefaultSetup(2);
screenNumber = max(Screen('Screens')); % use 2nd monitor if present
white = WhiteIndex(screenNumber); % color of fixation/response crosses
black = BlackIndex(screenNumber); % color of screen
cfix = ntsc2rgb([0.8 0 0]); % color of fixation target
cres = ntsc2rgb([0.8 .3 -.3]); % color of response target

%% LabJack setup
% Creates labjack object using labjack.m. Set silentMode to true in object
% construction ('silentMode',true) to run without communicating with LabJack or saving TTL/eye
% tracker data. Automatically enters silentMode if no LabJack connected.
daq = labJack(); %'verbose',true will include labjack.m text output

% LabJack FIO port settings
daqinfo.pulsePort = 0; %TTL pulse port output
daqinfo.xPort = 5; %x coordinate voltage input
daqinfo.yPort = 7; %y coordinate voltage input
daqinfo.diamPort = 3; %pupil diameter voltage input
daqinfo.pulsewidth = 10; %milliseconds
daqinfo.recBuff = zeros(16,1); % Read bytes (length = 10+ 2*Num of channels to read AIN)

% ConfigIO byte. See Labjack U3 documentation (5.2.3)
% Sets all channels other than FIO0 and FIO1 to analog input, Use ioBuffer[11] to change channel
% config (8 bit command, 0=digital,1=analog ex: 252 = 11111100)
if ~daq.silentMode
    ioBuffer = [0 248 3 11 0 0 13 0 64 0 252 0]';
    ioBuffer = daq.checksum(ioBuffer,'extended');
    bytesWritten = daq.rawWrite(ioBuffer);
    if bytesWritten == 0, warning('ConfigIO error: No bytes written.');end
    if bytesWritten < 12, warning('ConfigIO error: Did not write all of the buffer.');end
    recBuffer = zeros(12,1);
    [bytesRead,retBytes] = daq.rawRead(recBuffer,12);
    if bytesRead == 0, warning('ConfigIO error: Read failed.');end
    if bytesRead < 12, warning('ConfigIO error: Did not read all of the buffer.');end
end

% Bytes to send to read AIN from eye tracker channels
sendBuff = [0 248 5 0 0 0 0 1 daqinfo.xPort 31 1 daqinfo.yPort 31 1 daqinfo.diamPort 31]';
sendBuff = daq.checksum(sendBuff,'extended');

%% Position info
fixDim = 16;        % fixation target diameter (pixels)
picDim = 267;       % inter-target length (pixels)
textSize = 80;
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

%% Timing info
switch lower(delayParadigm)
    case 'fixed'
        fixRng = [1.0];
        iscontinuous = false;
    case 'discrete'
        fixRng = [0.5 1.5];
        iscontinuous = false;
    case 'continuous'
        fixRng = [0.5 1.5];
        isContinuous = true;
end
rspDurC = 1;        % response interval duration - CCDT (sec)
feedbackDur = 10;    % feedback display duration (sec)
stimTestInt = [0.750 1.250];    % Pulse interval during "Press any Key" screen
iti = [2.5 3.5];          % intertrial interval duration (sec)
ifi = Screen('GetFlipInterval', window); % inter-frame interval (s)
rspFramesC = round(rspDurC/ifi);
feedFrames = round(feedbackDur/ifi);
waitframes = 1; % number of frames to wait before re-drawing

%% Build Session Configuration
%sessConfig = struct();
sessConfig.cfix = cfix;
sessConfig.cres = cres;
sessConfig.fixDim = fixDim;
sessConfig.textSize = textSize;
sessConfig.rspDurC = rspDurC*1000;
%% Final setup
Ntpr = Ntarg; % Ntpr = targets per repetition
cont = 1; % Continue? boolean. Must start as true.
totalRT = 0;
[gpInd, gpName] = GetGamepadIndices;

Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);
vbl = Screen('Flip', window);

%% Test task/pulse loop. Send jittered pulses (750-1250ms) unitl keypress
keypress = false;
while ~keypress
    stimFrames = round((stimTestInt(1)+diff(stimTestInt)*rand)/ifi);
    for frame=1:stimFrames
        Screen('TextSize',window,60);
        DrawFormattedText(window,sprintf('Test screen.\n\nPlease wait for Administrator to begin.'),'center','center',white);
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        [keypress,~,keyCode] = KbCheck;
        if keyCode(KbName('ESCAPE')), sca; return; end
        if keypress, break;end
    end
    sendPulse;
end
HideCursor();

%% experimental loop
while cont
    for irep = 1:Ncal
        targ = randperm(Ntarg);
        itarg = 1;
        reactTimeMS = []; %initialize reation time array
        targetPix = [targetPix;[Px(targ),Py(targ)]];
        readInput(99);% Set target to 99 to indicate loop start
        
        while itarg<=min(Ntpr,Ntarg)
            respbit = 0; % reset response bit
            readInput(98);% Set target to 98 to indicate trial start
            
            %% intertrial interval
            itiDur = iti(1) + diff(iti)*rand;
            itiFrames = round(itiDur/ifi);
            for frame = 1:itiFrames
                if frame == 11 || frame == 21 || frame == 31, sendPulse; end % Must delay iti pulse to work with alignment
                readInput(0);
                vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                [~,~,keyCode] = KbCheck;
                if keyCode(KbName('ESCAPE')), respbit = 1; sca; return; end
            end
            
            %% fixation interval
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
                if frame == 1, trialStart = vbl; sendPulse; end% get stim onset in ms
                if frame== 11, sendPulse; end % 2 pulses on frame 1 and 11
                % Read response bit
                [~,~,keyCode] = KbCheck;
                if Gamepad('GetButton',gpInd,3) || keyCode(KbName('DownArrow'))
                    rspTime = GetSecs;
                    respbit = 1;
                    stimMS = (trialStart + fixDur)*1000;
                    trialRT = rspTime*1000-stimMS;
                    trialPressTime = rspTime*1000;
                    trialDelayTime = fixDur*1000;
                    break; % Restart trial on early press
                end
                if keyCode(KbName('ESCAPE')), respbit = 1; sca; return; end
            end
            if respbit
                reactTimeMS = [reactTimeMS;[trialRT,trialPressTime,trialDelayTime]];
                continue;
            end
            
            %% response interval
            for frame = 1:rspFramesC
                readInput(itarg+0.1);
                if frame == 1, sendPulse; end
                Screen('FillRect', window, cres, fixCoord(:,targ(itarg)));
                vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                % Read response bit
                [~,rspTime,keyCode] = KbCheck;
                if Gamepad('GetButton',gpInd,3) || keyCode(KbName('DownArrow'))
                    rspTime = GetSecs;
                    respbit = 1;
                    stimMS = (trialStart + fixDur)*1000;
                    trialRT = rspTime*1000-stimMS;
                    trialPressTime = rspTime*1000;
                    trialDelayTime = fixDur*1000;
                    totalRT = totalRT + trialRT;
                    break;
                end
                if keyCode(KbName('ESCAPE')), respbit = 1; sca; return; end
            end
            if respbit
                reactTimeMS = [reactTimeMS;[trialRT,trialPressTime,trialDelayTime]];
                itarg = itarg + 1;
            else
                reactTimeMS = [reactTimeMS;[rspDurC*1000,GetSecs*1000,fixDur*1000]]; %Save non-response
            end
        end
        %% Save data
        if exist(sessdir,'dir') ~=7 %check if first entry in session
            mkdir(sessdir);
        end
        cd (sessdir);
        if irep == 1, save('sessConfig.mat','sessConfig'); end
        sessRTs = [sessRTs;reactTimeMS];
        save('sessRTs.mat','sessRTs');
        save('syncTimes.mat','syncTimes');
        save('targetPix.mat','targetPix');
        save('eyeData.mat','eyeData');
        
        %% Feedback
        cd(recordsDir);
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
        
        %Set feedback text colors
        textColors = zeros(4,4);
        textColors(:,1) = [round(mean(feedbackRTs)); records.avgRT; round(min(feedbackRTs)); records.bestRT];
        for j=1:size(textColors,1)
            if textColors(j,1)<=250
                textColors(j,[2 3 4]) = [0 255 0];
            elseif textColors(j,1)<=500
                textColors(j,[2 3 4]) = [255 255 0];
            else
                textColors(j,[2 3 4]) = [255 0 0];
            end
        end

        for fFrame = 1:feedFrames
            %Header
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
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        end
    end
    %% End of set. Ask to continue. Send jittered pulses (750-1250ms)
    cont = 0;
    while ~cont
        stimFrames = round((stimTestInt(1)+diff(stimTestInt)*rand)/ifi);
        for frame=1:stimFrames
            Screen('TextSize',window,textSize);
            DrawFormattedText(window,'Press Button 3 to continue','center',screenYpixels*0.4,white);
            %Screen('TextSize',window,60);
            %DrawFormattedText(window,['esc to exit' '\nButton 3 to continue'],'center',screenYpixels*0.6,white);
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            [keypress,~,keyCode] = KbCheck;
            if Gamepad('GetButton',gpInd,3) || keyCode(KbName('DownArrow'))
                cont = 1;
            end
            if keyCode(KbName('ESCAPE')), sca; return; end
            if keypress, break;end
        end
        sendPulse;
    end
end
Priority(0);
sca;

%% Local functions
function sendPulse
% Send TTL pulse
if daq.silentMode, return; end
time = GetSecs*1000+(daqinfo.pulsewidth/2); %getting timestamp for pulse peak. Cmd-response for labjack is <1ms
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