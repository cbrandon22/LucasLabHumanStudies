function [gazeP, eeg,srate] = CCDTcalibrate(dataDir,session,liveDataFeed)
load(fullfile(dataDir,'behavioral',session,'eyeData'));
load(fullfile(dataDir,'behavioral',session,'targetPix'));
eeg = [];
srate=0;
if liveDataFeed
    % open jacksheet
    jacFile = fullfile(dataDir,'docs','jacksheet.txt');
    fid = fopen(jacFile,'r');
    JAC=textscan(fid,'%d%s','delimiter','\t');
    JAC{:,2} = strtrim(JAC{:,2});
    fclose(fid);
    channels = JAC{1};
    excludeChan = {'EKG1','EKG2','DC1','DC2','DC3','DC4','DC5','DC6','DC7'};
    [lia,locb] = ismember(excludeChan,JAC{2});
    if sum(lia)>0
        channels(locb(locb~=0)) = [];
    end
    flist = dir(fullfile(dataDir,'eeg.noreref'));
    sessChanBase = strsplit(flist(20).name,'.'); % assuming that the 20th file in this directory is a channel file
    sessChanBase = sessChanBase{1};
    fbase = fullfile(dataDir,'eeg.noreref',sessChanBase);
    c1eeg = look(fbase,channels(1),[],1)';
    eeg = ones(length(channels),size(c1eeg,2))*NaN;
    eeg(1,:) = c1eeg;
    for i=2:length(channels)
        eeg(i,:)=look(fbase,channels(i),[],1)';
    end
    fid = fopen(fullfile(dataDir,'eeg.noreref','params.txt'),'r');
    tline = fgetl(fid);
    while ischar(tline)
        splitline = strsplit(tline);
        if strcmp(splitline{1},'samplerate')
            srate = str2double(splitline{2});
        end
        tline=fgetl(fid);
    end
    fclose(fid);
end

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