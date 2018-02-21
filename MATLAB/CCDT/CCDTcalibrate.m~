function gazeP = CCDTcalibrate(subj,session)
dataDir = fullfile('~/Documents/CCDT/data',subj,session);
cd(dataDir);
load('eyeData');
load('targetPix');

% Calibration
% eyeData = [xVolt,yVolt,diamVolt,targetNumber](targ# 99 = trial start)
% targetPix = [xPixel, ypixel](ind=targetNumber)
% CC = new calibration constants
gazeP = zeros(length(eyeData),3);
startInd = find(eyeData(:,4)== 99); %Index of each trial start
endInd = [startInd-1;length(eyeData)];
endInd = endInd(endInd ~= 0); %Index of each trial end
N = size(unique(targetPix,'rows'),1); % Number of targets
Ncal = size(targetPix,1)/N;% Number of loops through each target
for i=1:Ncal % Calibrate voltages for each loop
    loopData = eyeData(startInd(i):endInd(i),[1 2 4]);
    V = loopData(loopData(:,3)>0 & loopData(:,3)< N+1,:); % XY Voltages and target during fixation
    mV = zeros(N,2);% Want a single x/y voltage pair for each target
    for ii = 1:N
        mV(ii,:) = median(V(V(:,3)==ii,1:2)); % robust Xvolt Yvolt for target
    end
    CC = [mV ones(N,1)]\targetPix(1+N*(i-1):N+N*(i-1),:); % least-squares solution of overdetermined linear system
    CC = CC'; % 2x3 matrix
    gazeP(startInd(i):endInd(i),:) = [eyeData(startInd(i):endInd(i),1)*CC(1,1) + eyeData(startInd(i):endInd(i),2)*CC(1,2) + CC(1,3),...
        eyeData(startInd(i):endInd(i),1)*CC(2,1) + eyeData(startInd(i):endInd(i),2)*CC(2,2) + CC(2,3), eyeData(startInd(i):endInd(i),4)];
end
end