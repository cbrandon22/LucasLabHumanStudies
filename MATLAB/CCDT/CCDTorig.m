function CCDT
% function CCDT
%   Color-change detection task

% color info
sca;
PsychDefaultSetup(2);
screenNumber = max(Screen('Screens')); % use 2nd monitor if present
white = WhiteIndex(screenNumber); % color of fixation/response crosses
black = BlackIndex(screenNumber); % color of screen
cfix = ntsc2rgb([0.8 0 0]); % color of fixation target
cres = ntsc2rgb([0.8 .3 -.3]); % color of response target

% position info
fixDim = 16;        % fixation target diameter (pixels)
picDim = 267;       % inter-target length (pixels)
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
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

% timing info
fixRng = [0.5 1.1]; % fixation interval range - CCDT (sec) %%% set instructed delay times here
rspDurC = 0.5;      % response interval duration - CCDT (sec)
iti = 0.5;          % intertrial interval duration (sec)
ifi = Screen('GetFlipInterval', window); % inter-frame interval (s)
rspFramesC = round(rspDurC/ifi);
itiFrames = round(iti/ifi);
waitframes = 1; % numer of frames to wait before re-drawing


% experimental loop
Ncal = 1; % Ncal = repetitions per calibration target
Ntpr = Ntarg; % Ntpr = targets per repetition
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);
vbl = Screen('Flip', window);
for irep = 1:Ncal
    targ = randperm(Ntarg);
    itarg = 1;
    while itarg<=min(Ntpr,Ntarg)
        respbit = 0; % reset response bit
        
        % intertrial interval
        for frame = 1:itiFrames
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        end
        
        % fixation interval
        fixDur = fixRng(1) + diff(fixRng)*rand;
        fixFrames = round(fixDur/ifi);
        for frame = 1:fixFrames
            Screen('FillRect', window, cfix, fixCoord(:,targ(itarg)));
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            % Read response bit
            [~,~,keyCode] = KbCheck;
            if keyCode(KbName('DownArrow')), respbit = 1; break; end
        end
        if respbit, continue; end
        
        % response interval
        for frame = 1:rspFramesC
            Screen('FillRect', window, cres, fixCoord(:,targ(itarg)));
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            % Read response bit
            [~,~,keyCode] = KbCheck;
            if keyCode(KbName('DownArrow')), respbit = 1; break; end
            if keyCode(KbName('ESCAPE')), respbit = 1; sca; return; end
        end
        if respbit, itarg = itarg + 1; end
    end
end

% end
Priority(0);
sca;