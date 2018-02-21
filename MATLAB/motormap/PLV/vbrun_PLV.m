function [plvAll, plvSourceSpace, plvBrodmannSpace, t] = vbrun_PLV(dat, ind, talStruct, resampleFreq, subj)

%% Description
% Calc triggered PLV between desired brain regions

% Inputs - dat: load from run_gete_ms (channels x trials x time)
%          ind: load from run_gete_ms
%          talStruct: loaded from anat_data for each subject (used to
%          create regionStruct)
%          resampleFreq: used for run_gete_ms (should be 500)
%          subj: cell array of subj name strings

% Outputs - plvAll: PLV pairwise across all channels
%           plvSourceSpace: PLV pairwise across all averaged data from
%           source space, defined by freesurfer pipieline (Ramayya et. al)
%           t: time axis for plotting PLV

% Workflow -- Run gete_ms, permute data to make chan x time x trials, then use 
%a binary fileMask to define move/ins/wait trials, 
% then run movePLV, insPLV, waitPLV. Then plot PLV between various channel
% combinations. Can also consider averaging PLV across that time interval
% for each electrode pair and then plotting an averaged move vs. ins. vs.
% wait PLV between each electrode pair. Then average the mean PLV across
% all source space pairs for wait vs. ins. vs. move.

% VPB - vivek.buch@uphs.upenn.edu
% University of Pennsylvania
% Dept of Neurosurgery, PGY3

%%
saveon=1;
ploton=0;
statson=0;
offset=-1; %in seconds, from eegdat file

allRegions=cell(1,length(talStruct));
allChans=zeros(1,length(talStruct));
% allBrodmann=cell(1,length(talStruct));

clear i
for i=1:length(allRegions)
    allRegions{i}=strrep(talStruct(i).avgSurf.anatRegion, ' ', '_'); %while creating variable, replace all spaces and hyphens with underscores, and add prefix for L vs. R
    allRegions{i}=strrep(allRegions{i}, '-', '_');
    allBrodmann{i}=strrep(talStruct(i).Loc5, ' ', '_'); %while creating variable, replace all spaces and hyphens with underscores, and add prefix for L vs. R
    allBrodmann{i}=strrep(allBrodmann{i}, '-', '_');
    allChans(i)=talStruct(i).channel;
    if istrue(talStruct(i).x<0);
        allRegions{i}=['L_' allRegions{i}];
        allBrodmann{i}=['L_' allBrodmann{i}];
    else
        allRegions{i}=['R_' allRegions{i}];
        allBrodmann{i}=['R_' allBrodmann{i}];
    end
    plvAll.name{i}=allRegions{i}; % for saving purposes
    plvAll.baname{i}=allBrodmann{i};
    plvAll.chans(i)=allChans(i);
end

unqReg=unique(allRegions);
unqBA=unique(allBrodmann);

clear i j currMask
j=0;
for i=1:length(unqReg);
    if length(unqReg{i}) <= 3;
        disp([unqReg{i} ' is not a brain region, skipping.'])
        continue
    end
    j=j+1;
    eval(sprintf('currMask=strcmp(''%s'',allRegions);',unqReg{i}));
    regionStruct.name{j}=unqReg{i};
    regionStruct.chans{j}=allChans(currMask);
    
    % for saving purposes
    plvSourceSpace.name{j}=unqReg{i};
    plvSourceSpace.chans{j}=allChans(currMask);
end

clear i j currMask
j=0;
for i=1:length(unqBA);
    if length(unqBA{i}) <= 3;
        disp(unqBA{i})
        continue
    end
    j=j+1;
    eval(sprintf('currMask=strcmp(''%s'',allBrodmann);',unqBA{i}));
    brodmannStruct.name{j}=unqBA{i};
    brodmannStruct.chans{j}=allChans(currMask);
    
    % for saving purposes
    plvBrodmannSpace.name{j}=unqBA{i};
    plvBrodmannSpace.chans{j}=allChans(currMask);
end

dat=permute(dat, [1 3 2]);
meanDat=zeros(length(regionStruct.chans),size(dat,2),size(dat,3));
brodmannDat=zeros(length(brodmannStruct.chans),size(dat,2),size(dat,3));

clear xy
for xy=1:size(meanDat,1);
    temp=mean(dat(regionStruct.chans{xy},:,:),1);
    meanDat(xy,:,:)=temp;
    clear temp
end

clear xy
for xy=1:size(brodmannDat,1);
    temp=mean(dat(brodmannStruct.chans{xy},:,:),1);
    brodmannDat(xy,:,:)=temp;
    clear temp
end

t=(0:size(meanDat,2)-1)/resampleFreq;
t=t+offset; %account for offset here


bands={'theta' 'beta' 'lowgamma' 'highgamma'};

for ij= 1:length(bands)
    freq=bands{ij};
    switch freq
        case 'theta'
            filtSpec.range=[3 12];
            filtSpec.order=250;
        case 'beta'
            filtSpec.range=[14 30];
            filtSpec.order=100;
        case 'lowgamma'
            filtSpec.range=[35 50];
            filtSpec.order=50;
        case 'highgamma'
            filtSpec.range=[70 100];
            filtSpec.order=50;
    end
    % filtspec.order depends on freq band --> 1000/(resampleFreq/filtSpec.order) >= 4*(1000/approx desired freq (Hz))
    disp(freq);
    % start running PLV
    % all individual channels pairwise
    disp('Running All-channel Pairwise PLV');
    
    eval(sprintf('plvAll.moveAll.%s=PLV(dat(:,:,ind.move), resampleFreq, filtSpec);',freq));
    eval(sprintf('plvAll.instructAll.%s=PLV(dat(:,:,ind.ins), resampleFreq, filtSpec);',freq));
    eval(sprintf('plvAll.waitAll.%s=PLV(dat(:,:,ind.wait), resampleFreq, filtSpec);',freq));
    
%     eval(sprintf('plvAll.moveL.%s=PLV(dat(:,:,ind.move&ind.left), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvAll.moveR.%s=PLV(dat(:,:,ind.move&ind.right), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvAll.instructL.%s=PLV(dat(:,:,ind.ins&ind.left), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvAll.isntructR.%s=PLV(dat(:,:,ind.ins&ind.right), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvAll.waitL.%s=PLV(dat(:,:,ind.wait&ind.left), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvAll.waitR.%s=PLV(dat(:,:,ind.wait&ind.right), resampleFreq, filtSpec);',freq));
%     
    % source-space averaged channels
    disp('Running Source Space PLV')
    
    eval(sprintf('plvSourceSpace.moveAll.%s=PLV(meanDat(:,:,ind.move), resampleFreq, filtSpec);',freq));
    eval(sprintf('plvSourceSpace.instructAll.%s=PLV(meanDat(:,:,ind.ins), resampleFreq, filtSpec);',freq));
    eval(sprintf('plvSourceSpace.waitAll.%s=PLV(meanDat(:,:,ind.wait), resampleFreq, filtSpec);',freq));
    
%     eval(sprintf('plvSourceSpace.moveL.%s=PLV(meanDat(:,:,ind.move&ind.left), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvSourceSpace.moveR.%s=PLV(meanDat(:,:,ind.move&ind.right), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvSourceSpace.instructL.%s=PLV(meanDat(:,:,ind.ins&ind.left), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvSourceSpace.isntructR.%s=PLV(meanDat(:,:,ind.ins&ind.right), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvSourceSpace.waitL.%s=PLV(meanDat(:,:,ind.wait&ind.left), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvSourceSpace.waitR.%s=PLV(meanDat(:,:,ind.wait&ind.right), resampleFreq, filtSpec);',freq));
    
    % brodmann area averaged channels
    disp('Running Brodmann Space PLV')
    
    eval(sprintf('plvBrodmannSpace.moveAll.%s=PLV(brodmannDat(:,:,ind.move), resampleFreq, filtSpec);',freq));
    eval(sprintf('plvBrodmannSpace.instructAll.%s=PLV(brodmannDat(:,:,ind.ins), resampleFreq, filtSpec);',freq));
    eval(sprintf('plvBrodmannSpace.waitAll.%s=PLV(brodmannDat(:,:,ind.wait), resampleFreq, filtSpec);',freq));
    
%     eval(sprintf('plvBrodmannSpace.moveL.%s=PLV(dat(:,:,ind.move&ind.left), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvBrodmannSpace.moveR.%s=PLV(dat(:,:,ind.move&ind.right), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvBrodmannSpace.instructL.%s=PLV(brodmannDat(:,:,ind.ins&ind.left), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvBrodmannSpace.isntructR.%s=PLV(brodmannDat(:,:,ind.ins&ind.right), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvBrodmannSpace.waitL.%s=PLV(brodmannDat(:,:,ind.wait&ind.left), resampleFreq, filtSpec);',freq));
%     eval(sprintf('plvBrodmannSpace.waitR.%s=PLV(brodmannDat(:,:,ind.wait&ind.right), resampleFreq, filtSpec);',freq));
    
end

%% plot;

if ploton
    figure;
    
    
    ROI{1}='L_Precentral_Gyrus';
    ROI{2}='L_Middle_Frontal_Gyrus';
    ROIind=find(strcmpi(ROI{1}, plvSourceSpace.name));
    ROIind(2)=find(strcmpi(ROI{2}, plvSourceSpace.name))
    
    ROI{1}='L_Brodmann_area_4_6';
    ROI{2}='L_Brodmann_area_1_2_3_5';
    baROIind=find(strcmpi(ROI{1}, plvBrodmannSpace.name));
    baROIind(2)=find(strcmpi(ROI{2}, plvBrodmannSpace.name))
    ROIind=find(strcmpi(ROI{2}, plvAll.baname))
    
    
    subplot(2,2,4)
    plot(t, squeeze(sourcespaceMovePLV_theta(:,1,2)), 'b');
    hold on
    plot(t, squeeze(sourcespaceIndPLV_theta(:,1,2)), 'r');
    plot(t, squeeze(sourcespaceWaitPLV_theta(:,1,2)), 'g');
    
    % title('Theta')
    
    legend move instruct wait
    xlabel('Post-stimulus time (s)')
    ylabel('PLV')
    title('L MST <--> MC connectivity: High gamma ')
    
    
    
    %mean figures with SEM 95% CI
    mu=mean(plvSourceSpace.move.lowgamma(:,7,12));
    mu(2)=mean(plvSourceSpace.instruct.lowgamma(:,7,12));
    
    sem=SEM_calc(plvSourceSpace.move.lowgamma(:,7,12), 0.05);
    sem(2)=SEM_calc(plvSourceSpace.instruct.lowgamma(:,7,12), 0.05);
    
    figure;
    H=bar(mu);
    set(H,'EdgeColor','b','FaceColor',[0.5,0.5,1],'LineWidth',1.5)
    set(gca,'XTickLabel',{'move'; 'instruct'})
    ylabel('mean PLV')
    title('Mean L PMC <--> MC connectivity: Low gamma ')
    xlim([0.5 length(mu)+0.5])
    
    hold on
    clear ii
    for ii=1:length(mu)
        plot([ii,ii],[mu(ii)-sem(ii),mu(ii)+sem(ii)],'-k','LineWidth',8)
    end
    hold off
    
    
end

%% save

if saveon
    save(['/media/vpbuch/VIVEK/motormap/' subj '_PLVnew'], 'regionStruct', 'brodmannStruct', 'plv*', 't');
end

%% Stats

if statson
    [tmp,p,Ci]=ttest2(plvSourceSpace.move.lowgamma(:,7,12), plvSourceSpace.instruct.lowgamma(:,7,12));
    
end
