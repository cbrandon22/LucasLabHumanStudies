%% Description

% Batch process gete_ms for multiple channels and subjects

% Inputs -  channels: cell array corresponding to each subj num channels
%           subj: cell array of each subject's name preceding '_events.mat' 
%           others same as gete_ms
% Output - dat: (channels x trials x time)
%          ind: trial mask for ind.move, ind.wait, ind.instruct

% VPB - vivek.buch@uphs.upenn.edu
% University of Pennsylvania
% Dept of Neurosurgery, PGY3


%% Define inputs
durationMS=3000;
offsetMS=-1000;
bufferMS=1000;
resampleFreq=512;
dir = '/media/vpbuch/VIVEK/motormap';

load('/media/vpbuch/VIVEK/motormap/allSubjNo88.mat')
subj=allSubj;

dbstop if error
%%

for j=1:length(subj)
    load([dir '/events/' subj{j} '_events.mat'])
    load([dir '/data/' subj{j} '/' subj{j} '_anatStruct.mat']) 
    disp(subj{j})
    
    type=cell(1, length(events));
%     item=cell(1, length(events));
    for ii=1:length(events);
        type{ii}=events(ii).type;
%         item{ii}=events(ii).item;
    end
    clear ii
    ind.move=strcmp('MOVE', type);
    ind.ins=strcmp('INSTRUCT', type);
    ind.wait=strcmp('WAIT', type);
%     
%     ind.left=strcmpi('Left Hand', item);
%     ind.right=strcmpi('Right Hand', item);
    
    
    ijwait=0;
    ijmove=0;
    for i=1:length(events)
        if istrue(ind.wait(i))
            ijwait=ijwait+1;
            for xx=1:length(anatStruct_mono)
                [temp] = gete_ms(anatStruct_mono(xx).elecNum, events(i), durationMS, offsetMS, bufferMS, [], [], [], resampleFreq);
%                 plot(temp)
%                 drawnow
                datWait(xx,ijwait,:)=temp;
                clear temp
            end
        elseif istrue(ind.move(i))
            ijmove=ijmove+1;
            for xy=1:length(anatStruct_mono)
                [temp] = gete_ms(anatStruct_mono(xy).elecNum, events(i), durationMS, offsetMS, bufferMS, [], [], [], resampleFreq);
%                 plot(temp)
%                 drawnow
                datMove(xy,ijmove,:)=temp;
                clear temp
            end
        else 
            continue
        end
    end
    clear i xx xy ijmove ijwait
    

    
    savedir='/media/vpbuch/VIVEK/motormap';
    save([savedir '/' subj{j} '_eegdat.mat'], 'dat*', 'ind') 
    clear dat* ind type
    disp('done')

end