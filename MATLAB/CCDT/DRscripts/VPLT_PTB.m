function VPLT_PTB
%   Visual preferential learning task
%   Uses Psychtoolbox-3 for experimental control and LabJack U3-LV for
%   inputs/outputs with TDT.
% 
%   DR 07/2016

% initialize LabJack
NET.addAssembly('LJUDDotNet'); % make the UD .NET assembly visible in MATLAB
ljObj = LabJack.LabJackUD.LJUD;
[~, ljH] = ljObj.OpenLabJack(LabJack.LabJackUD.DEVICE.U3, LabJack.LabJackUD.CONNECTION.USB, '0', true, 0); % open LabJack
ljObj.ePut(ljH, LabJack.LabJackUD.IO.PIN_CONFIGURATION_RESET, 0, 0, 0); % reset to factory defaults
ljObj.ePut(ljH, LabJack.LabJackUD.IO.PUT_ANALOG_ENABLE_PORT, 0, bin2dec('0000000000001111'), int32(16)); % configure flexible IO as analog (channel 0-3) and digital (channels 4-15) 

% picture set
debug = 0; % plot photodetector box? (0 or 1)
stimset = 50; % set number
stimdir = 'C:\TDT\Matlab\VPLT\stimsets\'; % set directory
cd([stimdir 'SET' sprintf('%03d',stimset)]); % open set directory
Npic = 200; % assumes 200 images
C = imread('1.bmp'); Cblack = zeros(size(C),'uint8'); Cpict = zeros([size(C),Npic],'uint8');
for ii = 1:Npic
    C = imread([num2str(ii) '.bmp']);
    Cpict(:,:,:,ii) = reshape(C,[size(C),1]);
end

% color info
sca;
PsychDefaultSetup(2);
screenNumber = max(Screen('Screens'));  % use 2nd monitor if present
white = WhiteIndex(screenNumber);       % color of fixation/response crosses
black = BlackIndex(screenNumber);       % color of screen
cfix = ntsc2rgb([0.8 0 0]);             % color of fixation target
cres = ntsc2rgb([0.8 .3 -.3]);          % color of response target

% position info
fixDim = 16;        % fixation target diameter (pixels; ~0.3 deg VA)
picDim = 267;       % half-width/height of pictures (pixels; ~11 deg VA)
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
[xCenter, yCenter] = RectCenter(windowRect);
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
pdtCoord = [20 screenYpixels-80 80 screenYpixels-20];
[Px,Py] = meshgrid(-1:1,-1:1); Ntarg = 9; % 9-point calibration
Px = Px(:)*picDim+xCenter; Py = Py(:)*picDim+yCenter;
baseRect = [0 0 fixDim fixDim];
fixCoord = NaN*ones(4,Ntarg);
for ii = 1:Ntarg
    fixCoord(:,ii) = CenterRectOnPointd(baseRect, Px(ii), Py(ii));
end
picRect = [0 0 2*picDim 2*picDim];
picCoord = CenterRectOnPointd(picRect, screenXpixels/2, screenYpixels/2);
fixCrossCoord = [-2*fixDim 2*fixDim 0 0; 0 0 -2*fixDim 2*fixDim];
fixCrossWidth = 5;

% timing info
fixRng = [0.5 1.1]; % fixation interval range - CCDT (sec)
fixDurV = 1.0;      % fixation interval duration - VPLT (sec)
rspDurC = 2.0;      % response interval duration - CCDT (sec)
rspDurV = 3.0;      % response interval duration - VPLT (sec)
iti = 1;          % intertrial interval duration (sec)
ifi = Screen('GetFlipInterval', window); % inter-frame interval (s)
fixFramesV = round(fixDurV/ifi);
rspFramesC = round(rspDurC/ifi);
rspFramesV = round(rspDurV/ifi);
itiFrames = round(iti/ifi);
waitframes = 1; % numer of frames to wait before re-drawing

% eye calibration
rep = 5; % Ntarget repetitions per calibration
CAL = NaN*ones(1000,10); ical = 1; % initialize calibration data matrix: rewarded trial x [Xpixel Ypixel Xvolt Yvolt CC(1,:) CC(2,:)], assumes no more than 1000 rewards per session
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);
vbl = Screen('Flip', window);
endbit = CCDT(rep,Ntarg);
if endbit, cleanup; return; end
CC = calibration(CAL,ical,Ntarg*rep);

% experimental loop
setsizes = [2 3 4]; % number of pictures in a set (each repeated twice)[3 4 5]
ipic = 1;
vbl = Screen('Flip', window);
while ipic<Npic
    setsize = min(setsizes(randi(3)),Npic-ipic+1);
    pics = repmat(ipic:ipic+setsize-1,1,2);
    io = randperm(length(pics));
    ii = 1;
    while ii<=length(pics)
        respbit = 0; % reset response bit
        
        % fixation interval
        ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.PUT_DIGITAL_PORT, 4, 201, 8, 0); ljObj.GoOne(ljH); % output word = 201
        ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.GET_AIN, 2, 0, 0, 0); % input eyeX
        ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.GET_AIN, 3, 0, 0, 0); % input eyeY
        for frame = 1:fixFramesV
            Screen('DrawLines', window, fixCrossCoord, fixCrossWidth, white, [xCenter yCenter], 2);
            if debug, Screen('FillRect', window, white, pdtCoord); end
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            ljObj.GoOne(ljH);
            [~,eyeX] = ljObj.GetResult(ljH, LabJack.LabJackUD.IO.GET_AIN, 2, 0);
            [~,eyeY] = ljObj.GetResult(ljH, LabJack.LabJackUD.IO.GET_AIN, 3, 0);
            Xp = eyeX*CC(1,1)+eyeY*CC(1,2)+CC(1,3); % convert current gaze location from volts to pixels
            Yp = eyeX*CC(2,1)+eyeY*CC(2,2)+CC(2,3);
            if Xp<picCoord(1) || Xp>picCoord(3) || Yp<picCoord(2) || Yp>picCoord(4) % must only fixate within picCoord -- could narrow this
                respbit = 1; break; 
            end
        end
        if respbit, continue; end % stay in fixation interval until fixates for fixDurV sec
        
        % response interval
        ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.PUT_DIGITAL_PORT, 4, pics(io(ii)), 8, 0); ljObj.GoOne(ljH); % output word = 1-200
        ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.GET_AIN, 2, 0, 0, 0); % input eyeX
        ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.GET_AIN, 3, 0, 0, 0); % input eyeY
        cpic = squeeze(Cpict(:,:,:,pics(io(ii))));
        picTexture = Screen('MakeTexture', window, cpic);
        for frame = 1:rspFramesV
            Screen('DrawTexture', window, picTexture, [], picCoord);
            if debug, Screen('FillRect', window, white, pdtCoord); end
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            ljObj.GoOne(ljH);
            [~,eyeX] = ljObj.GetResult(ljH, LabJack.LabJackUD.IO.GET_AIN, 2, 0);
            [~,eyeY] = ljObj.GetResult(ljH, LabJack.LabJackUD.IO.GET_AIN, 3, 0);
            Xp = eyeX*CC(1,1)+eyeY*CC(1,2)+CC(1,3); % convert current gaze location from volts to pixels
            Yp = eyeX*CC(2,1)+eyeY*CC(2,2)+CC(2,3);
            if Xp<picCoord(1) || Xp>picCoord(3) || Yp<picCoord(2) || Yp>picCoord(4)
                respbit = 1; break; % remove picture if gaze not on it
            end
        end
        ii = ii + 1;
        
        % intertrial interval
        ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.PUT_DIGITAL_PORT, 4, 0, 8, 0); ljObj.GoOne(ljH); % output word = 0
        ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.GET_DIGITAL_BIT, 13, 0, 0, 0); % input end bit (check only during intertrial)
        for frame = 1:itiFrames
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            ljObj.GoOne(ljH);
            [~,endbit] = ljObj.GetResult(ljH, LabJack.LabJackUD.IO.GET_DIGITAL_BIT, 13, 0);
            if endbit, cleanup; return; end
        end
    end
    disp(pics(io)); % display the set of pics for this loop
    ipic = ipic + setsize;
    
    % CCDT
    endbit = CCDT(1,5);
    if endbit, cleanup; return; end
    CC = calibration(CAL,ical,Ntarg*rep);
end
cleanup;

    % nested color-change detection task (CCDT) function 
    function endbit = CCDT(Ncal,Ntpr) % Ncal = repetitions per calibration target, Ntpr = targets per repetition
        endbit = 0;
        for irep = 1:Ncal
            targ = randperm(Ntarg);
            itarg = 1;
            while itarg<=min(Ntpr,Ntarg)
                respbit = 0; % reset response bit
                
                % intertrial interval
                ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.PUT_DIGITAL_PORT, 4, 210, 8, 0); ljObj.GoOne(ljH); % output word = 210
                ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.GET_DIGITAL_BIT, 13, 0, 0, 0); % input end bit (check only during intertrial)
                for frame = 1:itiFrames
                    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                    ljObj.GoOne(ljH);
                    [~,endbit] = ljObj.GetResult(ljH, LabJack.LabJackUD.IO.GET_DIGITAL_BIT, 13, 0);
                    if endbit, return; end
                end
                
                % fixation interval
                ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.PUT_DIGITAL_PORT, 4, targ(itarg)+210, 8, 0); ljObj.GoOne(ljH); % output word = 211-219
                ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.GET_DIGITAL_BIT, 12, 0, 0, 0); % input joystick bit
                fixDur = fixRng(1) + diff(fixRng)*rand;
                fixFrames = round(fixDur/ifi);
                for frame = 1:fixFrames
                    Screen('FillRect', window, cfix, fixCoord(:,targ(itarg)));
                    if debug, Screen('FillRect', window, white, pdtCoord); end
                    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                    ljObj.GoOne(ljH);
                    [~,respbit] = ljObj.GetResult(ljH, LabJack.LabJackUD.IO.GET_DIGITAL_BIT, 12, 0);
                    if respbit, break; end
                end
                if respbit, continue; end
                
                % response interval
                ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.PUT_DIGITAL_PORT, 4, 220, 8, 0); ljObj.GoOne(ljH); % output word = 220
                ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.GET_DIGITAL_BIT, 12, 0, 0, 0); % input joystick bit
                for frame = 1:rspFramesC
                    Screen('FillRect', window, cres, fixCoord(:,targ(itarg)));
                    if debug, Screen('FillRect', window, white, pdtCoord); end
                    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                    ljObj.GoOne(ljH);
                    [~,respbit] = ljObj.GetResult(ljH, LabJack.LabJackUD.IO.GET_DIGITAL_BIT, 12, 0);
                    if respbit, break; end
                end
                if respbit
                    ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.GET_AIN, 2, 0, 0, 0); % input eyeX
                    ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.GET_AIN, 3, 0, 0, 0); % input eyeY
                    ljObj.GoOne(ljH);
                    [~,eyeX] = ljObj.GetResult(ljH, LabJack.LabJackUD.IO.GET_AIN, 2, 0);
                    [~,eyeY] = ljObj.GetResult(ljH, LabJack.LabJackUD.IO.GET_AIN, 3, 0);
                    CAL(ical,:) = [Px(targ(itarg)), Py(targ(itarg)), eyeX, eyeY, CC(1,:), CC(2,:)]; % calibration info for each rewarded CCDT trial
                    ical = ical + 1;
                    itarg = itarg + 1;
                end
            end
        end
        
        % intertrial interval
        ljObj.AddRequest(ljH, LabJack.LabJackUD.IO.PUT_DIGITAL_PORT, 4, 210, 8, 0); ljObj.GoOne(ljH); % output word = 210
        for frame = 1:itiFrames
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        end
    end

    % nested cleanup function
    function cleanup
        Priority(0);
        sca;
        save(['C:\TDT\Matlab\VPLT\cal_' datestr(now,'mmddyy-HHMM') '.mat'],'CAL','picCoord','-mat');
    end
end

function CC = calibration(CAL,ical,ncal)
% Calibration function
% CAL = [Xpixel Ypixel Xvolt Yvolt CC(1,:) CC(2,:)]
% ical = last+1 data into CAL
% ncal = number of data points to use for calibration
% CC = new calibration constants
ind = max(ical-ncal,1):(ical-1);
P = CAL(ind,1:2);
V = CAL(ind,3:4);
[uP,~,iu] = unique(P,'rows');
N = size(uP,1);
mV = zeros(N,2);
for ii = 1:N
    mV(ii,:) = median(V(iu==ii,:)); % robust Xvolt Yvolt for target
end
CC = [mV ones(N,1)]\uP; % least-squares solution of overdetermined linear system
CC = CC'; % 2x3 matrix
end