%% stplv plots from humanCCDT_connectomic_st saved output
clear
% Input
allSubj = {'HUP069','HUP133','HUP136','HUP139','HUP140',...
    'HUP142','HUP143','HUP145','HUP146','HUP150'};%,'HUP152','HUP153',...
    %'HUP154','HUP157'};
plvdir = '/Volumes/HumanStudies/HumanStudies/CCDT/scratch/connectivity';

% Loop through each subject
for s=1:length(allSubj)
    load(fullfile(plvdir,[allSubj{s} '_Wst0.mat']));
    disp(allSubj{s})
    disp(length(W_st(1).precue))
%     for i=1:length(vRT)
%         for f=1:length(W_st(1).precue)
%             % sample graph metrics calculation (strength and clustering coef)
%             str(i,:,f) = strengths_und(W_st(i).precue(f).adj);
%             cCa(i,:,f) = clustering_coef_wu(W_st(i).precue(f).adj);
%         end
%     end
%     keyboard
%     %sample scatter for trial-by-trial RT as function of strength
%     figure
%     scatter(nanmean(str(vRT>0&vDT>1450,:,1),2),vRT(vRT>0&vDT>1450));
%     figure
%     scatter(nanmean(str(vRT>0&vDT>1450,:,2),2),vRT(vRT>0&vDT>1450));
%     
%     %sample plot for average strength for fastest trials
%     figure
%     plot(mean(str(ifast500,:,2),1));
%     hold on
%     plot(mean(str(islow500,:,2),1));
%     
%     figure
%     plot(mean(str(ilate500&vRT>0,:,2),1));
%     hold on
%     plot(mean(str(iearly500&vRT>0,:,2),1));
end