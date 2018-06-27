%% set directories and set ieeg portal session name
pwPath = '/Volumes/HumanStudies/MATLAB/LucasLabHumanStudies/MATLAB/cbr_ieeglogin.bin';
stimDir = '/Volumes/HumanStudies/HumanStudies/corticalStim';
portalSessName = 'HUP128_phaseII_D02';

%% Initialize session, search for subject stim info
sess = IEEGSession(portalSessName,'cbrandon22',pwPath);
subj = strsplit(sess.data.snapName,'_');
subj = subj{1};
load(fullfile(stimDir,'stim_masterSubjInfo.mat'));
sInd = length(stim_masterSubjInfo);

if ~ismember(subj,{stim_masterSubjInfo.subj})
    disp('Use ieeg portal to add this subject session info')
    sInd = sInd+1;
    keyboard
    stim_masterSubjInfo(sInd).subj = subj;
    stim_masterSubjInfo(sInd).sessStart = 440222;
    stim_masterSubjInfo(sInd).sessEnd = 445077;
    save(fullfile(stimDir,'stim_masterSubjInfo.mat'),'stim_masterSubjInfo');
end

cLbls = sess.data.channelLabels(:,1);
excludeChan = {'EKG1','EKG2'};
inclCh = sess.data.rawChannels(~ismember(cLbls,excludeChan));
dataset = sess.data.deriveDataset([sess.userName '_' subj '_stim'],inclCh);