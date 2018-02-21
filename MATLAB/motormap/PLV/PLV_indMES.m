%Calculate hedges' g for all primary and association cortex, saving
%[hedgesg CI(1) CI(2)]in plvDir/individual
plvDir = '/Volumes/LUCAS_DRIVE/HumanStudies/motormap/scratch/figs/PLV/mes/anat/hedgesg/lowFreq250_highFreq50';
load(fullfile(plvDir,'allMotor.mat'))
iterations = 3000; %number of samples for bootstrap CI
trimPLV = 512; %ms to trim from beginning and end of PLV
sampleWin = 512; %window size (samples)
sampleStep = 256; %window step size (samples)
cueInd = 1525-trimPLV; %Index of wait or go cue

m_hfMoveWait = zeros(size(m_hfMoveDat,1),3);
m_hfMoveMove = zeros(size(m_hfMoveDat,1),3);
m_hfWaitWait = zeros(size(m_hfMoveDat,1),3);
m_lfMoveWait = zeros(size(m_lfMoveDat,1),3);
m_lfMoveMove = zeros(size(m_lfMoveDat,1),3);
m_lfWaitWait = zeros(size(m_lfMoveDat,1),3);

for i=1:size(m_hfMoveDat,1)
    % High frequency motor
    moveDat = m_hfMoveDat(i,:);
    waitDat = m_hfWaitDat(i,:);
    moveDat = moveDat(trimPLV:(length(moveDat)-2*trimPLV))';
    waitDat = waitDat(trimPLV:(length(waitDat)-2*trimPLV))';
    tit = ['hf motor pair: ' num2str(i)];
    [moveWait,moveMove,waitWait] = PLV_mm_mesPlot(moveDat,waitDat,t,sampleWin,sampleStep,cueInd,tit,iterations);
    m_hfMoveWait(i,:) = moveWait;
    m_hfMoveMove(i,:) = moveMove;
    m_hfWaitWait(i,:) = waitWait;
    close;
    % Low frequency motor
    moveDat = m_lfMoveDat(i,:);
    waitDat = m_lfWaitDat(i,:);
    moveDat = moveDat(trimPLV:(length(moveDat)-2*trimPLV))';
    waitDat = waitDat(trimPLV:(length(waitDat)-2*trimPLV))';
    tit = ['lf motor pair: ' num2str(i)];
    [moveWait,moveMove,waitWait] = PLV_mm_mesPlot(moveDat,waitDat,t,sampleWin,sampleStep,cueInd,tit,iterations);
    m_lfMoveWait(i,:) = moveWait;
    m_lfMoveMove(i,:) = moveMove;
    m_lfWaitWait(i,:) = waitWait;
    close;
end
mkdir(fullfile(plvDir,'individual'))
cd(fullfile(plvDir,'individual'))
save('hfMotorStats', 'm_hfMoveWait','m_hfMoveMove','m_hfWaitWait')
save('lfMotorStats', 'm_lfMoveWait','m_lfMoveMove','m_lfWaitWait')
