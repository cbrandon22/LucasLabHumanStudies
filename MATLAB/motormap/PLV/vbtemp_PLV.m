% load eegdat and anatStruct for the patient
 
%  'theta'
filtSpecT.range=[3 12];
filtSpecT.order=250;
%  'beta'
filtSpecB.range=[14 30];
filtSpecB.order=100;
% 'lowgamma'
filtSpecLG.range=[35 50];
filtSpecLG.order=50;
% 'highgamma'
filtSpecHG.range=[70 200];
filtSpecHG.order=50;
 
datMove=permute(datMove, [1 3 2]);
datWait=permute(datWait, [1 3 2]);
 
plvAllMove_highgamma=PLV(datMove,512,filtSpecHG);
plvAllWait_highgamma=PLV(datWait,512,filtSpecHG);
 
plvAllMove_lowgamma=PLV(datMove,512,filtSpecLG);
plvAllWait_lowgamma=PLV(datWait,512,filtSpecLG);
 
plvAllMove_beta=PLV(datMove,512,filtSpecB);
plvAllWait_beta=PLV(datWait,512,filtSpecB);
 
plvAllMove_theta=PLV(datMove,512,filtSpecT);
plvAllWait_theta=PLV(datWait,512,filtSpecT);
 
t=(0:size(datMove,2)-1)/512;
t=t-1; %account for offset in seconds
 
%%
 
for i=1:length(anatStruct_mono); 
    regions{i}=anatStruct_mono(i).anatAbbr; 
end
 
ROI{1}='L-PRE';
ROI{2}='L-POST';
ROI{3}='L-F2';
ROI{4}='L-SMA';
ROI{5}='R-PRE';
ROI{6}='R-POST';
ROI{7}='R-F2';
ROI{8}='R-SMA';
 
ROI{9}='L-T3';
ROI{10}='L-HIP';
 
 
ROIind1=find(strcmpi(ROI{1}, regions));
ROIind2=find(strcmpi(ROI{2}, regions));
ROIind3=find(strcmpi(ROI{3}, regions));
ROIind4=find(strcmpi(ROI{4}, regions));
ROIind5=find(strcmpi(ROI{5}, regions));
ROIind6=find(strcmpi(ROI{6}, regions));
ROIind7=find(strcmpi(ROI{7}, regions));
ROIind8=find(strcmpi(ROI{8}, regions));
 
ROIind9=find(strcmpi(ROI{9}, regions));
ROIind10=find(strcmpi(ROI{10}, regions));
 
 
%% Now load in plvdat and run this per subject

currROIindOne=ROIind1;
currROIindTwo=ROIind2;
 
clear i j
count=0;
for i=1:length(currROIindOne)
    for j=1:length(currROIindTwo)
        count=count+1;
        fprintf('%d <--> %d ; ', currROIindOne(i), currROIindTwo(j));
        if currROIindOne(i)<currROIindTwo(j)
            plvMove_highgamma_POST_PRE(count,:)=plvAllMove_highgamma(:,currROIindOne(i),currROIindTwo(j));
            plvWait_highgamma_POST_PRE(count,:)=plvAllWait_highgamma(:,currROIindOne(i),currROIindTwo(j));
        else
            plvMove_highgamma_POST_PRE(count,:)=plvAllMove_highgamma(:,currROIindTwo(j),currROIindOne(i));
            plvWait_highgamma_POST_PRE(count,:)=plvAllWait_highgamma(:,currROIindTwo(j),currROIindOne(i));
        end
    end
end
clear i j count

%
% move epochs
avgPlvMove_highgamma_POST_PRE=mean(plvMove_highgamma_POST_PRE,1);
stdPlvMove_highgamma_POST_PRE=std(plvMove_highgamma_POST_PRE);
 
avgBaselineMove_highgamma_POST_PRE=mean(plvMove_highgamma_POST_PRE(:,250:500)); % average across channels
avgBaselineMove_highgamma_POST_PRE=mean(avgBaselineMove_highgamma_POST_PRE); % average entire baseline period
stdBaselineMove_highgamma_POST_PRE=std(mean(plvMove_highgamma_POST_PRE(:,250:500),2));
 
 
normPlvMove_highgamma_POST_PRE=(avgPlvMove_highgamma_POST_PRE-avgBaselineMove_highgamma_POST_PRE)./stdBaselineMove_highgamma_POST_PRE;
 
 
% wait epochs
avgPlvWait_highgamma_POST_PRE=mean(plvWait_highgamma_POST_PRE,1);
stdPlvWait_highgamma_POST_PRE=std(plvWait_highgamma_POST_PRE);
 
avgBaselineWait_highgamma_POST_PRE=mean(plvWait_highgamma_POST_PRE(:,250:500)); % average across channels
avgBaselineWait_highgamma_POST_PRE=mean(avgBaselineWait_highgamma_POST_PRE); % average entire baseline period
stdBaselineWait_highgamma_POST_PRE=std(mean(plvWait_highgamma_POST_PRE(:,250:500),2));
 
 
normPlvWait_highgamma_POST_PRE=(avgPlvWait_highgamma_POST_PRE-avgBaselineWait_highgamma_POST_PRE)./stdBaselineWait_highgamma_POST_PRE;
 
%% plots
 
% figure;
% errorbar(t, fgsmooth(normPlvMove_highgamma_F2_PRE,10), stdPlvMove_highgamma_F2_PRE, 'r')
% hold on
% plot(t, fgsmooth(normPlvMove_highgamma_F2_PRE,10), 'r', 'linewidth', 2)
% xlim([-0.5 2])
