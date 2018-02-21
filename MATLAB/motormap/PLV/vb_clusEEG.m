
%% Cluster-based EEG

% Workflow
% 1. Assign each electrode to a particular cluster (from Cameron) for each
% condition (moveWait, moveInstruct, instructWait)
% 2. Get the ECOG signal for each electrode and
% create a matrix of ECOG signal channel x time for all MOVE trials, WAIT
% trials, INSTRUCT trials within each cluster. 


%% Define inputs (from defineClusters.m)

% subj='HUP111'; 
% clus1chan=[];
% clus2chan=[];
% clus3chan=[];
% allClus={'clus1' 'clus2' 'clus3'};


%%
task = 'motormap';
dirs = le_dirs(task);
evDir = fullfile(dirs.events);
cd(evDir);
load([subj, '_events.mat']);
if useBipolar
    anatStruct_bipol = le_centroids2Anatstruct(subj);
else
    [~,anatStruct_mono] = le_centroids2Anatstruct(subj);
end

% load(['/media/vpbuch/VIVEK/motormap/events/' subj '_events.mat']);
% if useBipolar
%     load(['/media/vpbuch/VIVEK/motormap/data/' subj '/' subj '_anatStruct_bipol.mat'])
% else
%     load(['/media/vpbuch/VIVEK/motormap/data/' subj '/' subj '_anatStruct_mono.mat'])
% end

durationMS=10000;
offsetMS=-3000;
bufferMS=1000;
resampleFreq=512;
baseSecBetEv = 10;
baseJitt = 5;
basePriorMS = 500;
basePostMS = 500;

type=cell(1, length(events));
for ii=1:length(events);
type{ii}=events(ii).type;
end
clear ii
ind.move=strcmp('MOVE', type);
ind.ins=strcmp('INSTRUCT', type);
ind.wait=strcmp('WAIT', type);
baseEv=get_baseline_random(events,baseSecBetEv,baseJitt);

if useBipolar
    elecVec=zeros(length(anatStruct_bipol),2);
    for jj=1:length(elecVec)
        elecVec(jj,:)=anatStruct_bipol(jj).elecNum;
    end
    clear jj   
else   %using monopolar 
    clus1chan=unique(clus1chan);
    clus2chan=unique(clus2chan);
    clus0chan=unique(clus0chan);
    if exist('clusHchan','var')
        clusHchan=unique(clusHchan);
    end
    
    elecVec=zeros(1, length(anatStruct_mono));
    for jj=1:length(elecVec)
        elecVec(jj)=anatStruct_mono(jj).elecNum;
    end
    clear jj
end

alreadystarted=who('clus1Struct');
if isempty(alreadystarted)
    disp('preallocate clusters')
    clus1Struct=struct;
    clus2Struct=struct;
    clus0Struct=struct;
    clusHStruct=struct;
end

% add to preexisting clusters if they exist (actual cluster, not just
% preallocated above)
if isfield(clus1Struct, 'subj')
    disp('adding to existing clusters')
    stackval1=length(clus1Struct);
    stackval2=length(clus2Struct);
    stackval0=length(clus0Struct);
    if isfield(clusHStruct, 'subj')
        stackvalH=length(clusHStruct);
    else
        stackvalH=0;
    end
else
    stackval1=0;
    stackval2=0;
    stackval0=0;
    stackvalH=0;
end

disp([subj ' cluster-specific EEG...'])
if ~useBipolar %monopolar
    for j=1:length(allClus);
        currClus=allClus{j};
        switch currClus
            
            case 'clus1'
                chan=clus1chan;
                disp(currClus)
                stackval=stackval1;
                for i = 1:length(chan)
                    disp(chan(i))
                    [eegch] = gete_ms(chan(i), events, durationMS, offsetMS, bufferMS, [], [], [], resampleFreq);
                    [baseEEG] = gete_ms(chan(i), baseEv, basePriorMS+basePostMS, -basePriorMS, bufferMS, [], [], [], resampleFreq);
                    clus1Struct(i+stackval).subj=subj;
                    clus1Struct(i+stackval).elec=chan(i);
                    clus1Struct(i+stackval).datMove=eegch(ind.move,:);
                    clus1Struct(i+stackval).datWait=eegch(ind.wait,:);
                    clus1Struct(i+stackval).datBase=baseEEG;
                    clus1Struct(i+stackval).elecLbl=anatStruct_mono(elecVec==chan(i)).elecLbl;
                    clus1Struct(i+stackval).anatLbl=anatStruct_mono(elecVec==chan(i)).anatAbbr;
                    clear eegch
                end
                clear i
                
            case 'clus2'
                chan=clus2chan;
                disp(currClus)
                stackval=stackval2;
                for i = 1:length(chan)
                    disp(chan(i))
                    [eegch] = gete_ms(chan(i), events, durationMS, offsetMS, bufferMS, [], [], [], resampleFreq);
                    [baseEEG] = gete_ms(chan(i), baseEv, basePriorMS+basePostMS, -basePriorMS, bufferMS, [], [], [], resampleFreq);
                    clus2Struct(i+stackval).subj=subj;
                    clus2Struct(i+stackval).elec=chan(i);
                    clus2Struct(i+stackval).datMove=eegch(ind.move,:);
                    clus2Struct(i+stackval).datWait=eegch(ind.wait,:);
                    clus2Struct(i+stackval).datBase=baseEEG;
                    clus2Struct(i+stackval).elecLbl=anatStruct_mono(elecVec==chan(i)).elecLbl;
                    clus2Struct(i+stackval).anatLbl=anatStruct_mono(elecVec==chan(i)).anatAbbr;
                    clear eegch
                end
                clear i
                
            case 'clus0'
                chan=clus0chan;
                disp(currClus)
                stackval=stackval0;
                for i = 1:length(chan)
                    disp(chan(i))
                    [eegch] = gete_ms(chan(i), events, durationMS, offsetMS, bufferMS, [], [], [], resampleFreq);
                    [baseEEG] = gete_ms(chan(i), baseEv, basePriorMS+basePostMS, -basePriorMS, bufferMS, [], [], [], resampleFreq);
                    clus0Struct(i+stackval).subj=subj;
                    clus0Struct(i+stackval).elec=chan(i);
                    clus0Struct(i+stackval).datMove=eegch(ind.move,:);
                    clus0Struct(i+stackval).datWait=eegch(ind.wait,:);
                    clus0Struct(i+stackval).datBase=baseEEG;
                    clus0Struct(i+stackval).elecLbl=anatStruct_mono(elecVec==chan(i)).elecLbl;
                    clus0Struct(i+stackval).anatLbl=anatStruct_mono(elecVec==chan(i)).anatAbbr;
                    clear eegch
                end
                clear i
                
            case 'clusH'
                if exist('clusHchan','var')
                    chan=clusHchan;
                    disp(currClus)
                    stackval=stackvalH;
                    for i = 1:length(chan)
                        disp(chan(i))
                        [eegch] = gete_ms(chan(i), events, durationMS, offsetMS, bufferMS, [], [], [], resampleFreq);
                        [baseEEG] = gete_ms(chan(i), baseEv, basePriorMS+basePostMS, -basePriorMS, bufferMS, [], [], [], resampleFreq);
                        clusHStruct(i+stackval).subj=subj;
                        clusHStruct(i+stackval).elec=chan(i);
                        clusHStruct(i+stackval).datMove=eegch(ind.move,:);
                        clusHStruct(i+stackval).datWait=eegch(ind.wait,:);
                        clusHStruct(i+stackval).datBase=baseEEG;
                        clusHStruct(i+stackval).elecLbl=anatStruct_mono(elecVec==chan(i)).elecLbl;
                        clusHStruct(i+stackval).anatLbl=anatStruct_mono(elecVec==chan(i)).anatAbbr;
                        clear eegch
                    end
                end
                clear i

                
        end
    end
    clear j
end


clearvars -except clus1Struct clus2Struct clus0Struct clusHStruct useBipolar clus_info subjVec allClus

disp('done')


