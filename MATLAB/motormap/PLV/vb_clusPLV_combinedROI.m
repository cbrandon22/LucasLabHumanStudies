%% Cluster-based PLV
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

%% Define inputs (load subj_clusPLV_PW_clusx)
subjList = {'HUP084','HUP087','HUP091','HUP092','HUP111'};
%subjList = {'HUP091'};
for sub=1:length(subjList)
    for jj=1:length(allClus)
        freqLabels = {'lowFreq250_highFreq50'};
        for thisFreq=1:length(freqLabels)
            subj=subjList{sub};
            currClus=allClus{jj};
            dirs = le_dirs('motormap');
            plvDir = fullfile(dirs.scratch,'PHASE/k2/btwRegions',freqLabels{thisFreq});
            saveDir = fullfile(dirs.scratch,'PHASE/k2/btwRegions',freqLabels{thisFreq});
            cd(plvDir);
            disp(['Looking for: ' subj ': ' freqLabels{thisFreq} ' ' currClus])
            if exist([subj '_clusPLV_PW_' currClus '.mat'],'file')==2
                load([subj '_clusPLV_PW_' currClus]);
                disp(['loaded: ' subj ': ' freqLabels{thisFreq} ' ' currClus])
                
                %%
                motorROI={'L-PRE' 'L-POST' 'R-PRE' 'R-POST'};
                motorROIind = find(ismember(regions,motorROI));
                nonMotorROIind = find(~ismember(regions,motorROI));
                % nonMotorROI will be defined below as everything else
                % hippROI={'L-HIP' 'L-PHIP' 'L-AMYG' 'R-HIP' 'R-PHIP' 'R-AMYG'};
                
                % motorROIind=NaN(length(motorROI),length(regions));
                % hippROIind=NaN(length(hippROI),length(regions));
                
%                 clear i
%                 for i=1:length(motorROI)
%                     temp=find(strcmpi(motorROI(i),regions));
%                     if ~isempty(temp)
%                         motorROIind(i,1:length(temp))=temp;
%                     end
%                     clear temp
%                 end
%                 motorROIind=unique(motorROIind(~isnan(motorROIind)))';
%
%                 nonMotorROIind =(1:1:length(regions));
%                 clear i
%                 for i=1:length(motorROIind)
%                     nonMotorROIind(motorROIind(i))=NaN;
%                 end
%                 nonMotorROIind=unique(nonMotorROIind(~isnan(nonMotorROIind)));
                
                
                
                % clear i
                % for i=1:length(hippROI)
                %     temp=find(strcmpi(hippROI(i),regions));
                %     if ~isempty(temp)
                %         hippROIind(i,1:length(temp))=temp;
                %     end
                %     clear temp
                % end
                % hippROIind=unique(hippROIind(~isnan(hippROIind)))';
                
                %%
                
                % all pairwise motor and non-motor
                clear i j
                count=0;
                for i=motorROIind
                    for j=motorROIind
                        if i<j
                            count=count+1;
                            motor_normThetaMove(count,:)=m2_normThetaMove(:,i,j)';
                            motor_normThetaWait(count,:)=m2_normThetaWait(:,i,j)';
                            motor_normGammaMove(count,:)=m2_normGammaMove(:,i,j)';
                            motor_normGammaWait(count,:)=m2_normGammaWait(:,i,j)';
                            motor_ThetaMove(count,:)=plvThetaMove(:,i,j)';
                            motor_ThetaWait(count,:)=plvThetaWait(:,i,j)';
                            motor_GammaMove(count,:)=plvGammaMove(:,i,j)';
                            motor_GammaWait(count,:)=plvGammaWait(:,i,j)';
                        else
                            continue
                        end
                    end
                end
                
                clear i j
                count=0;
                for i=nonMotorROIind
                    for j=nonMotorROIind
                        if i<j
                            count=count+1;
                            nonMotor_normThetaMove(count,:)=m2_normThetaMove(:,i,j)';
                            nonMotor_normThetaWait(count,:)=m2_normThetaWait(:,i,j)';
                            nonMotor_normGammaMove(count,:)=m2_normGammaMove(:,i,j)';
                            nonMotor_normGammaWait(count,:)=m2_normGammaWait(:,i,j)';
                            nonMotor_ThetaMove(count,:)=plvThetaMove(:,i,j)';
                            nonMotor_ThetaWait(count,:)=plvThetaWait(:,i,j)';
                            nonMotor_GammaMove(count,:)=plvGammaMove(:,i,j)';
                            nonMotor_GammaWait(count,:)=plvGammaWait(:,i,j)';
                        else
                            continue
                        end
                    end
                end
                btwROI = [motorROIind, nonMotorROIind];
                clear i j
                count=0;
                for i=nonMotorROIind
                    for j=motorROIind
                        if i<j
                            count=count+1;
                            btw_normThetaMove(count,:)=m2_normThetaMove(:,i,j)';
                            btw_normThetaWait(count,:)=m2_normThetaWait(:,i,j)';
                            btw_normGammaMove(count,:)=m2_normGammaMove(:,i,j)';
                            btw_normGammaWait(count,:)=m2_normGammaWait(:,i,j)';
                            btw_ThetaMove(count,:)=plvThetaMove(:,i,j)';
                            btw_ThetaWait(count,:)=plvThetaWait(:,i,j)';
                            btw_GammaMove(count,:)=plvGammaMove(:,i,j)';
                            btw_GammaWait(count,:)=plvGammaWait(:,i,j)';
                        else
                            continue
                        end
                    end
                end
                clear i j count
                
                % clear i j
                % count=0;
                % for i=hippROIind
                %     for j=hippROIind
                %         if i<j
                %             count=count+1;
                %             hipp_normThetaMove(count,:)=m2_normThetaMove(:,i,j)';
                %             hipp_normThetaWait(count,:)=m2_normThetaWait(:,i,j)';
                %             hipp_normGammaMove(count,:)=m2_normGammaMove(:,i,j)';
                %             hipp_normGammaWait(count,:)=m2_normGammaWait(:,i,j)';
                %         else
                %             continue
                %         end
                %     end
                % end
                % clear i j count
                
                %%
                % save
                saveon=1;
                if saveon
                    disp(['saving: ' subj '_' currClus ' in ' saveDir])
                    cd(saveDir)
                    save([subj sprintf('_clusPLV_combinedROI_%s.mat',currClus)], 'motor*', 'nonMotor*', 'btw*', 'regions','t') % 'hipp*'
                end
            end
        end
        clearvars -except clus1Struct clus2Struct clus0Struct clusHStruct useBipolar clus_info subjVec allClus subjList sub
        disp('done')
    end
end