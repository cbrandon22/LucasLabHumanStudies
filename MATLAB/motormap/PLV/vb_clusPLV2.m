%% Cluster-based PLV
% *****Added random baseline rather than pre-cue (difference from vb_clusPLV)*****
% 1. What is the connectivity between the middle frontal,
% pre-central, post-central, and temporal electrodes WITHIN each cluster?
% Then compare this pattern between clusters.
% 2. Connectivity between those region electrodes BETWEEN clusters?
% 3. Connectivity within a single region from electrodes belonging to the
% same vs. different clusters?

%% Workflow
% 1. Assign each electrode to a particular cluster (from Cameron) for each
% condition (moveWait, moveInstruct, instructWait)
% 2. Get the ECOG signal for each electrode and
% create a matrix of ECOG signal channel x time for all MOVE trials, WAIT
% trials, INSTRUCT trials within each cluster. 

% Steps 1 & 2 done in vb_clusEEG.m

% 3. Compute PLV for each pair of electrodes within each cluster per
% subject so can be normalized and compared across patients
% 4. Normalize each PLV to its own baseline period
% 5a. Average the normalized PLV for all electrode pairs across subjects within and between
% each brain region (designated as "m2_" below)

% 5b. Vs. Average raw PLV within brain regions, then normalize (designated
% "plv..." below) -- this method does not allow for averaging all PLV
% signals across subjects, drops it down to one normalized PLV per condition per
% subject

%% Define inputs
subjList = {'HUP084','HUP087','HUP091','HUP092','HUP111'};
for sub=1:length(subjList)
    for jj=1:2
        freqStruct.labels = {'broadTheta250_broadTheta300';'alpha200_beta100';'delta600_theta250';'lowFreq250_highFreq50'};
        freqStruct.range1 = [3,12;8,12;1,4;3,33];
        freqStruct.range2 = [3,12;12,32;4,8;70,100];
        freqStruct.order1 = [250;200;600;250];
        freqStruct.order2 = [300;100;250;50];
        for thisFreq=1:length(freqStruct.labels)
            subj=subjList{sub};
            currClus=jj;
            dirs = le_dirs;
            % loadVar='media/vpbuch/VIVEK/motormap/clusterEEGdat_clus0_Long_allSubjMono_moveWait.mat'; % confirm condition, created from vb_clusEEG.m
            saveon=1;
            saveDir = fullfile(dirs.scratch,'PHASE/k2/btwRegions',freqStruct.labels{thisFreq});
            useBipolar=0;
            offset=-3; %seconds
            baseSecBetEv = 10;
            baseJitt = 5;
            basePriorMS = 500;
            basePostMS = 500;
            
            %% Sensor-space PLV analysis
            % create currStruct
            alreadyloaded=who('clus1Struct');
            if isempty(alreadyloaded)
                load(loadVar)
            end
            
            eval(sprintf('currClusStruct=clus%dStruct;',currClus));
            clusSubjs = {currClusStruct.subj};
            currStructInd = strcmp(subj,clusSubjs);
            currStruct = currClusStruct(currStructInd);
            
            clear i
            for i=1:length(currStruct)
                datMove(i,:,:)=currStruct(i).datMove;
                datWait(i,:,:)=currStruct(i).datWait;
                datBase(i,:,:)=currStruct(i).datBase;
            end
            currDatMove=permute(datMove, [1 3 2]);
            currDatWait=permute(datWait, [1 3 2]);
            currDatBase=permute(datBase, [1 3 2]);
            
            filtSpecT.range=freqStruct.range1(thisFreq,:);
            filtSpecT.order=freqStruct.order1(thisFreq);
            filtSpecG.range=freqStruct.range2(thisFreq,:);
            filtSpecG.order=freqStruct.order2(thisFreq);
            
            disp([subj ': running cluster-based PLV analysis... cluster ' num2str(currClus) ', freqs: ' freqStruct.labels{thisFreq}])
            
            [plvThetaMove]=PLV(currDatMove, 512, filtSpecT);
            [plvThetaWait]=PLV(currDatWait, 512, filtSpecT);
            [plvThetaBase]=PLV(currDatBase, 512, filtSpecT);
            [plvGammaMove]=PLV(currDatMove, 512, filtSpecG);
            [plvGammaWait]=PLV(currDatWait, 512, filtSpecG);
            [plvGammaBase]=PLV(currDatBase, 512, filtSpecG);
            
            %pre-allocate loop variable for speed
            m2_normThetaMove=zeros(size(currDatMove,2),size(currDatMove,1),size(currDatMove,1));
            m2_normThetaWait=zeros(size(currDatMove,2),size(currDatMove,1),size(currDatMove,1));
            m2_normGammaMove=zeros(size(currDatMove,2),size(currDatMove,1),size(currDatMove,1));
            m2_normGammaWait=zeros(size(currDatMove,2),size(currDatMove,1),size(currDatMove,1));
            
            
            clear i j
            for i=1:size(plvGammaMove,2);
                for j=1:size(plvGammaMove,3);
                    m2_normGammaMove(:,i,j)=((plvGammaMove(:,i,j))-mean(plvGammaBase(:,i,j),1))./std(plvGammaBase(:,i,j));
                    m2_normThetaMove(:,i,j)=((plvThetaMove(:,i,j))-mean(plvThetaBase(:,i,j),1))./std(plvThetaBase(:,i,j));
                    m2_normGammaWait(:,i,j)=((plvGammaWait(:,i,j))-mean(plvGammaBase(:,i,j),1))./std(plvGammaBase(:,i,j));
                    m2_normThetaWait(:,i,j)=((plvThetaWait(:,i,j))-mean(plvThetaBase(:,i,j),1))./std(plvThetaBase(:,i,j));
                end
            end
            
            % all pairwise
            clear i j
            count=0;
            for i=1:size(m2_normThetaMove,2)
                for j=1:size(m2_normThetaMove,3)
                    if i<j
                        count=count+1;
                        m2_PW_normThetaMove(count,:)=m2_normThetaMove(:,i,j)';
                        m2_PW_normThetaWait(count,:)=m2_normThetaWait(:,i,j)';
                        m2_PW_normGammaMove(count,:)=m2_normGammaMove(:,i,j)';
                        m2_PW_normGammaWait(count,:)=m2_normGammaWait(:,i,j)';
                        PW_ThetaMove(count,:)=plvThetaMove(:,i,j)';
                        PW_ThetaWait(count,:)=plvThetaWait(:,i,j)';
                        PW_GammaMove(count,:)=plvGammaMove(:,i,j)';
                        PW_GammaWait(count,:)=plvGammaWait(:,i,j)';
                    else
                        continue
                    end
                end
            end
            clear i j count
            
            m2_avg_PW_normThetaMove=mean(m2_PW_normThetaMove);
            m2_avg_PW_normThetaWait=mean(m2_PW_normThetaWait);
            m2_avg_PW_normGammaMove=mean(m2_PW_normGammaMove);
            m2_avg_PW_normGammaWait=mean(m2_PW_normGammaWait);
            
            
            
            %% ROI source-space analysis
            regions=cell(1,length(currStruct));
            clear i
            for i=1:length(currStruct)
                regions{i}=currStruct(i).anatLbl;
            end
            
            % ROI{1}='L-PRE';
            % ROI{2}='L-POST';
            % ROI{3}='L-F2';
            % ROI{4}='L-T3';
            % ROI{5}='L-HIP';
            % ROI{6}='L-PHIP';
            % ROI{7}='L-SMG';
            % ROI{8}='L-P2';
            % ROI{9}='R-PRE';
            % ROI{10}='R-POST';
            % ROI{11}='R-F2';
            % ROI{12}='R-T3';
            % ROI{13}='R-HIP';
            % ROI{14}='R-PHIP';
            % ROI{15}='R-SMG';
            % ROI{16}='R-P2';
            %
            % motorROI={'L-PRE' 'L-POST' 'L-F2' 'R-PRE' 'R-POST' 'R-F2'};
            % nonMotorROI={'L-HIP' 'L-PHIP' 'R-HIP' 'R-PHIP'};
            %
            % ROIind=NaN(length(ROI),length(currStruct));
            % motorROIind=NaN(length(motorROI),length(currStruct));
            % nonMotorROIind=NaN(length(nonMotorROI),length(currStruct));
            %
            % % find indices corresponding to above ROI, store a binary vector if its
            % % good vs. empty, and then rename the ROIs withot hyphens for saving purposes
            % clear i
            % for i=1:length(ROI)
            %     temp=find(strcmpi(ROI(i),regions));
            %     if ~isempty(temp)
            %         ROIind(i,1:length(temp))=temp;
            %     end
            %     goodROI(i)=length(unique(~isnan(ROIind(i,:))))-1; %report 1 if there are nonNan numbers, 0 otherwise
            %     ROI{i}=strrep(ROI{i},'-','');
            %     clear temp
            % end
            %
            % clear i
            % for i=1:length(motorROI)
            %     temp=find(strcmpi(motorROI(i),regions));
            %     if ~isempty(temp)
            %         motorROIind(i,1:length(temp))=temp;
            %     end
            %     clear temp
            % end
            % motorROIind=unique(motorROIind(~isnan(motorROIind)))';
            %
            % clear i
            % for i=1:length(nonMotorROI)
            %     temp=find(strcmpi(nonMotorROI(i),regions));
            %     if ~isempty(temp)
            %         nonMotorROIind(i,1:length(temp))=temp;
            %     end
            %     clear temp
            % end
            % nonMotorROIind=unique(nonMotorROIind(~isnan(nonMotorROIind)))';
            %
            %
            % %%
            % disp(ROI)
            % disp(goodROI)
            % inptOne=input('ROI one: ');
            % if isempty(inptOne)
            %     ROIon=0;
            %     ROIchk=0;
            % else
            %     currROIindOne=ROIind(inptOne,~isnan(ROIind(inptOne,:)));
            %     inptTwo=input('ROI two: ');
            %     currROIindTwo=ROIind(inptTwo,~isnan(ROIind(inptTwo,:)));
            %     ROIone=ROI{inptOne};
            %     ROItwo=ROI{inptTwo};
            %     ROIchk=1;
            % end
            %
            % %confirm that subject has channels in these regions, if not input another ROI
            % while ROIchk
            %     if isempty(currROIindOne)
            %         disp(goodROI)
            %         newOne=input('ROI one empty, please give new region number  \n (If desired regions empty, leave blank and hit Enter): __ ');
            %         if isempty(newOne) % desired region empty, then hit return
            %             disp('skipping ROI analysis')
            %             ROIon=0;
            %             ROIchk=0;
            %             continue
            %         else
            %             currROIindOne=ROIind(newOne,~isnan(ROIind(newOne,:)));
            %             ROIone=ROI{newOne}; %rename the string for saving purposes
            %             continue
            %         end
            %     elseif isempty(currROIindTwo)
            %         newTwo=input('ROI two empty, please give new region number  \n (If desired regions empty, leave blank and hit Enter): __');
            %         if isempty(newTwo)
            %             disp('skipping ROI analysis')
            %             ROIon=0;
            %             ROIchk=0;
            %             continue
            %         else
            %             currROIindTwo=ROIind(newTwo,~isnan(ROIind(newTwo,:)));
            %             ROItwo=ROI{newTwo}; % rename string
            %             continue
            %         end
            %     else
            %         ROIon=1;
            %         ROIchk=0;
            %     end
            % end
            %
            % if ROIon
            %     disp('running ROI source-space PLV analysis...')
            %     % average PLV across brain regions then normalize
            %
            %     %pre-allocate loop variables for speed
            %     plvThetaMove_ROI=zeros(length(currROIindOne)*length(currROIindTwo), size(currDatMove,2));
            %     plvThetaWait_ROI=zeros(length(currROIindOne)*length(currROIindTwo), size(currDatMove,2));
            %     plvGammaMove_ROI=zeros(length(currROIindOne)*length(currROIindTwo), size(currDatMove,2));
            %     plvGammaWait_ROI=zeros(length(currROIindOne)*length(currROIindTwo), size(currDatMove,2));
            %
            %     clear i j count
            %     count=0;
            %     for i=1:length(currROIindOne)
            %         for j=1:length(currROIindTwo)
            %             if currROIindOne(i)<currROIindTwo(j)
            %                 count=count+1;
            %                 fprintf('%d <--> %d ; ', currROIindOne(i), currROIindTwo(j));
            %                 plvThetaMove_ROI(count,:)=plvThetaMove(:,currROIindOne(i),currROIindTwo(j));
            %                 plvThetaWait_ROI(count,:)=plvThetaWait(:,currROIindOne(i),currROIindTwo(j));
            %                 plvGammaMove_ROI(count,:)=plvGammaMove(:,currROIindOne(i),currROIindTwo(j));
            %                 plvGammaWait_ROI(count,:)=plvGammaWait(:,currROIindOne(i),currROIindTwo(j));
            %             elseif currROIindOne(i)>currROIindTwo(j)
            %                 count=count+1;
            %                 fprintf('%d <--> %d ; ', currROIindOne(i), currROIindTwo(j));
            %                 plvThetaMove_ROI(count,:)=plvThetaMove(:,currROIindTwo(j),currROIindOne(i));
            %                 plvThetaWait_ROI(count,:)=plvThetaWait(:,currROIindTwo(j),currROIindOne(i));
            %                 plvGammaMove_ROI(count,:)=plvGammaMove(:,currROIindTwo(j),currROIindOne(i));
            %                 plvGammaWait_ROI(count,:)=plvGammaWait(:,currROIindTwo(j),currROIindOne(i));
            %             elseif currROIindOne(i)==currROIindTwo(j)
            %                 continue
            %             end
            %         end
            %     end
            %     clear i j count
            %
            %
            %     avg_plvThetaMove_ROI=mean(plvThetaMove_ROI,1);
            %     avgBase_plvThetaMove_ROI=mean2(plvThetaMove_ROI(:,baseTimeInd1:baseTimeInd2));
            %     stdBase_plvThetaMove_ROI=std2(plvThetaMove_ROI(:,baseTimeInd1:baseTimeInd2));
            %     %avgBase_plvThetaMove_ROI=mean(avg_plvThetaMove_ROI(baseTimeInd1:baseTimeInd2));
            %     %stdBase_plvThetaMove_ROI=std(avg_plvThetaMove_ROI(baseTimeInd1:baseTimeInd2));
            %     avg_plvThetaWait_ROI=mean(plvThetaWait_ROI,1);
            %     avgBase_plvThetaWait_ROI=mean2(plvThetaWait_ROI(:,baseTimeInd1:baseTimeInd2));
            %     stdBase_plvThetaWait_ROI=std2(plvThetaWait_ROI(:,baseTimeInd1:baseTimeInd2));
            %
            %     avg_plvGammaMove_ROI=mean(plvGammaMove_ROI,1);
            %     avgBase_plvGammaMove_ROI=mean2(plvGammaMove_ROI(:,baseTimeInd1:baseTimeInd2));
            %     stdBase_plvGammaMove_ROI=std2(plvGammaMove_ROI(:,baseTimeInd1:baseTimeInd2));
            %
            %     avg_plvGammaWait_ROI=mean(plvGammaWait_ROI,1);
            %     avgBase_plvGammaWait_ROI=mean2(plvGammaWait_ROI(:,baseTimeInd1:baseTimeInd2));
            %     stdBase_plvGammaWait_ROI=std2(plvGammaWait_ROI(:,baseTimeInd1:baseTimeInd2));
            %
            %     avg_normThetaMove_ROI=(avg_plvThetaMove_ROI-avgBase_plvThetaMove_ROI)./stdBase_plvThetaMove_ROI;
            %     avg_normThetaWait_ROI=(avg_plvThetaWait_ROI-avgBase_plvThetaWait_ROI)./stdBase_plvThetaWait_ROI;
            %     avg_normGammaMove_ROI=(avg_plvGammaMove_ROI-avgBase_plvGammaMove_ROI)./stdBase_plvGammaMove_ROI;
            %     avg_normGammaWait_ROI=(avg_plvGammaWait_ROI-avgBase_plvGammaWait_ROI)./stdBase_plvGammaWait_ROI;
            %
            %
            %
            %     % normalize each PLV then average across regions (just to compare to method 1 above)
            %
            %     %pre-allocate loop variables for speed
            %     m2_normThetaMove_ROI=zeros(length(currROIindOne)*length(currROIindTwo), size(currDatMove,2));
            %     m2_normThetaWait_ROI=zeros(length(currROIindOne)*length(currROIindTwo), size(currDatMove,2));
            %     m2_normGammaMove_ROI=zeros(length(currROIindOne)*length(currROIindTwo), size(currDatMove,2));
            %     m2_normGammaWait_ROI=zeros(length(currROIindOne)*length(currROIindTwo), size(currDatMove,2));
            %     clear i j count
            %     count=0;
            %     disp('running normalized ROI PLV analysis...')
            %     for i=1:length(currROIindOne)
            %         for j=1:length(currROIindTwo)
            %             if currROIindOne(i)<currROIindTwo(j)
            %                 count=count+1;
            %                 fprintf('%d <--> %d ; ', currROIindOne(i), currROIindTwo(j));
            %                 m2_normThetaMove_ROI(count,:)=m2_normThetaMove(:,currROIindOne(i),currROIindTwo(j));
            %                 m2_normThetaWait_ROI(count,:)=m2_normThetaWait(:,currROIindOne(i),currROIindTwo(j));
            %                 m2_normGammaMove_ROI(count,:)=m2_normGammaMove(:,currROIindOne(i),currROIindTwo(j));
            %                 m2_normGammaWait_ROI(count,:)=m2_normGammaWait(:,currROIindOne(i),currROIindTwo(j));
            %             elseif currROIindOne(i)>currROIindTwo(j)
            %                 count=count+1;
            %                 fprintf('%d <--> %d ; ', currROIindOne(i), currROIindTwo(j));
            %                 m2_normThetaMove_ROI(count,:)=m2_normThetaMove(:,currROIindTwo(j),currROIindOne(i));
            %                 m2_normThetaWait_ROI(count,:)=m2_normThetaWait(:,currROIindTwo(j),currROIindOne(i));
            %                 m2_normGammaMove_ROI(count,:)=m2_normGammaMove(:,currROIindTwo(j),currROIindOne(i));
            %                 m2_normGammaWait_ROI(count,:)=m2_normGammaWait(:,currROIindTwo(j),currROIindOne(i));
            %             elseif currROIindOne(i)==currROIindTwo(j)
            %                 continue
            %             end
            %         end
            %     end
            %     clear i j count
            %
            %     m2_avg_normThetaMove_ROI=mean(m2_normThetaMove_ROI,1);
            %     m2_avg_normThetaWait_ROI=mean(m2_normThetaWait_ROI,1);
            %     m2_avg_normGammaMove_ROI=mean(m2_normGammaMove_ROI,1);
            %     m2_avg_normGammaWait_ROI=mean(m2_normGammaWait_ROI,1);
            %
            % else
            %     disp('not enough desired ROIs for source-space analysis, skipping...')
            %     avg=0; %for saving purposes
            %
            % end
            %
            t=(0:size(currDatMove,2)-1)/512;
            t=t+offset;
            %
            %
            %
            %
            % % define save variable and save
            % saveVar=sprintf('%s_clus%d',subj,currClus);
            %
            % % eval(sprintf('%s_plvThetaMove=plvThetaMove;',saveVar));
            % % eval(sprintf('%s_plvThetaWait=plvThetaWait;',saveVar));
            % % eval(sprintf('%s_plvGammaMove=plvGammaMove;',saveVar));
            % % eval(sprintf('%s_plvGammaWait=plvGammaWait;',saveVar));
            % %
            % % eval(sprintf('%s_m2_normThetaMove=m2_normThetaMove;',saveVar));
            % % eval(sprintf('%s_m2_normThetaWait=m2_normThetaWait;',saveVar));
            % % eval(sprintf('%s_m2_normGammaMove=m2_normGammaMove;',saveVar));
            % % eval(sprintf('%s_m2_normGammaWait=m2_normGammaWait;',saveVar));
            %
            % if ROIon
            %     eval(sprintf('plvThetaMove_%s_%s=plvThetaMove_ROI;',ROIone,ROItwo));
            %     eval(sprintf('plvThetaWait_%s_%s=plvThetaWait_ROI;',ROIone,ROItwo));
            %     eval(sprintf('plvGammaMove_%s_%s=plvGammaMove_ROI;',ROIone,ROItwo));
            %     eval(sprintf('plvGammaWait_%s_%s=plvGammaWait_ROI;',ROIone,ROItwo));
            %
            %     eval(sprintf('avg_plvThetaMove_%s_%s=avg_plvThetaMove_ROI;',ROIone,ROItwo));
            %     eval(sprintf('avg_plvThetaWait_%s_%s=avg_plvThetaWait_ROI;',ROIone,ROItwo));
            %     eval(sprintf('avg_plvGammaMove_%s_%s=avg_plvGammaMove_ROI;',ROIone,ROItwo));
            %     eval(sprintf('avg_plvGammaWait_%s_%s=avg_plvGammaWait_ROI;',ROIone,ROItwo));
            %
            %     eval(sprintf('avg_normThetaMove_%s_%s=avg_normThetaMove_ROI;',ROIone,ROItwo));
            %     eval(sprintf('avg_normThetaWait_%s_%s=avg_normThetaWait_ROI;',ROIone,ROItwo));
            %     eval(sprintf('avg_normGammaMove_%s_%s=avg_normGammaMove_ROI;',ROIone,ROItwo));
            %     eval(sprintf('avg_normGammaWait_%s_%s=avg_normGammaWait_ROI;',ROIone,ROItwo));
            %
            %     eval(sprintf('m2_normThetaMove_%s_%s=m2_normThetaMove_ROI;',ROIone,ROItwo));
            %     eval(sprintf('m2_normThetaWait_%s_%s=m2_normThetaWait_ROI;',ROIone,ROItwo));
            %     eval(sprintf('m2_normGammaMove_%s_%s=m2_normGammaMove_ROI;',ROIone,ROItwo));
            %     eval(sprintf('m2_normGammaWait_%s_%s=m2_normGammaWait_ROI;',ROIone,ROItwo));
            %
            %     eval(sprintf('m2_avg_normThetaMove_%s_%s=m2_avg_normThetaMove_ROI;',ROIone,ROItwo));
            %     eval(sprintf('m2_avg_normThetaWait_%s_%s=m2_avg_normThetaWait_ROI;',ROIone,ROItwo));
            %     eval(sprintf('m2_avg_normGammaMove_%s_%s=m2_avg_normGammaMove_ROI;',ROIone,ROItwo));
            %     eval(sprintf('m2_avg_normGammaWait_%s_%s=m2_avg_normGammaWait_ROI;',ROIone,ROItwo));
            %
            % end
            %
            % disp('done')
            %%
            % save
            if saveon
                disp('saving...')
                if ~(exist(saveDir,'dir')==7),mkdir(saveDir);end
                cd(saveDir);
                save([subj sprintf('_clusPLV_PW_clus%d.mat',currClus)], 'm2*', 'plv*','PW*','regions','t')  %'avg*','ROIind*','goodROI', 'ROI'
            end
        end
            clearvars -except clus1Struct clus2Struct useBipolar clus_info subjVec allClus subjList sub
            disp('done')
    end
end