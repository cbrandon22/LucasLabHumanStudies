% CCDT Human Connectomic Analysis

%% Inputs
ddir = '/Volumes/HumanStudies/HumanStudies/CCDT'; % should contain events/ and eeg/ folders
saveDir = '/Volumes/HumanStudies/HumanStudies/CCDT/scratch/connectivity';
singleWstSubj = {'HUP069','HUP133','HUP136','HUP139','HUP140',...
     'HUP142','HUP143','HUP145','HUP146','HUP150'}; % single Wst0 file subjects
splitWstSubj = {'HUP152','HUP153','HUP154','HUP157','HUP160','HUP165','HUP168','HUP171'}; % subjects to split
% allSubj = {'HUP069','HUP133','HUP136','HUP139','HUP140',...
%      'HUP142','HUP143','HUP145','HUP146','HUP150','HUP152','HUP153','HUP154',...
%      'HUP160','HUP165','HUP168','HUP157','HUP171'};
allSubj = {'HUP160','HUP165','HUP168'};
saveon=1;
bipol=0; % 1 if want to use bipolar montage, 0 if monopolar
lpc=200; % low-pass cut-off frequency for butterworth filter: wn=lpc/(fs/2)
fs=512; % sample rate
twin=0.5; % seconds for single-trial PLV sliding window
saveDir = ['/Volumes/HumanStudies/HumanStudies/CCDT/scratch/connectivity/' num2str(twin*1000) 'msWin'];
%trial window settings
durationMS=3000;
offsetMS=-1000;
bufferMS=1000;
resampleFreq=512;
bT = .300; %(s) buffer time after end of DT

%filter specifications for PLV
filtSpec(1).range = [3 12];
filtSpec(1).order = 250;
filtSpec(1).name = 'theta';
filtSpec(2).range = [70 100];
filtSpec(2).order = 50;
filtSpec(2).name = 'high gamma';
filtSpec(3).range = [35 55];
filtSpec(3).order = 50;
filtSpec(3).name = 'low gamma';

%% CCDT Connectomics subject loop
for iji=1:length(allSubj)
    subj = allSubj{iji}; %subject name string vector
    if ismember(subj,splitWstSubj)
        splitTrials = 2; % split files by trial # (total trials/splitTrials = # of trials per W_st file). Only use when #channels>130
    else
        splitTrials = 1;
    end
    disp(subj)
    load(fullfile(ddir, 'events',[subj '_events.mat']));
    
    fid = fopen(fullfile(ddir,'eeg',subj, 'docs/jacksheet.txt'));
    C = textscan(fid,'%s');
    fclose(fid);
    JAC{:,1} = C{1}(1:2:end-1);
    JAC{:,2} = C{1}(2:2:end);

    chrec = cellfun(@str2double, JAC{1});
    chlbl = JAC{2};
    
    type=cell(1, length(events));
    vRT = NaN(length(events),1); % vector of RTs per trial over sessions
    vDT = NaN(length(events),1); % vector of DTs per trial over sessions
    for ii=1:length(events)
        type{ii}=events(ii).type;
        vDT(ii) = events(ii).delay;
        vRT(ii) = events(ii).rt;
    end
    clear ii
    ind.cue=strcmp('FIX_START', type);
    ind.go=strcmp('CC', type);
    ind.resp=strcmp('RESPONSE', type);
    
    D = struct;
    N = length(events(ind.cue));
    vDT = vDT(ind.cue);
    vRT = vRT(ind.cue);
    
    % data
    disp('Processing raw signal...')
    bpC = 0;
    for jj = 1:length(chrec)
        [dat] = gete_ms(chrec(jj), events, durationMS, offsetMS, bufferMS, [], [], [], resampleFreq);
        if bipol
            try %#ok<UNRCH>
                [datref] = gete_ms(chrec(jj+1), events, durationMS, offsetMS, bufferMS, [], [], [], resampleFreq);
                dat = dat - datref;
                bpC = bpC+1;
            catch
                disp('End of bipolar data')
                dat = [];
            end
        end
        
        if ~isempty(dat)
            [b,a] = butter(2,lpc/(fs/2));
            dat = filtfilt(b,a,dat);
            D(jj).raw = dat; %Raw dat per channel aligned on each event (cue, go, and response)
            D(jj).dat = dat(ind.cue,:);  % trials x time just aligned on cue
            %         for ix = 1:N
            %             D(jj).dat(ix,t>=(DT(ix)+(bT*1000))) = NaN; %throw out time after color change (+100ms buffer)
            %         end
            if ~bipol
                D(jj).ch = chrec(jj);
            else
                D(jj).ch = [num2str(chrec(jj)) '-' num2str(chrec(jj+1))];
            end
            
        end
        
        
    end
    
    t = ((1:size(D(1).dat,2))/resampleFreq)*1000;
    t = t+offsetMS;
    
    if ~bipol
        cDat = zeros(length(chrec),length(t), length(vDT)); %ch x time x trials
        for xix = 1:length(chrec)
            cDat(xix,:,:) = D(xix).dat';
        end
    else
        cDat = zeros(bpC,length(t), length(vDT)); %ch x time x trials
        for xix = 1:bpC
            cDat(xix,:,:) = D(xix).dat';
        end
    end
    
    if twin == 0.5
        tBind1 = find(t>=-1000,1)+1:find(t>-500,1)-1; %-1000 to -500ms baseline
        tBind = find(t>=-500,1)+1:find(t>=0,1)-1; %-500 to 0ms baseline
        tDTind1 = find(t>=0,1)+1:find(t>=500,1)-1; %0 to delay period
        tDTind2 = find(t>=500,1)+1:find(t>=1000,1)-1; %0 to delay period
        tDTind3 = find(t>=1000,1)+1:find(t>=1500,1)-1; %0 to delay period  
    elseif twin == 0.25
        tBind3 = find(t>=-1000,1)+1:find(t>-750,1)-1; %-1000 to -750ms baseline
        tBind2 = find(t>=-750,1)+1:find(t>=-500,1)-1; %-750 to -500ms baseline
        tBind1 = find(t>=-500,1)+1:find(t>=-250,1)-1; %-500 to -250ms baseline
        tBind = find(t>=-250,1)+1:find(t>=0,1)-1; %-250 to 0ms baseline
        tDTind1 = find(t>=0,1)+1:find(t>=250,1)-1; %0 to delay period
        tDTind2 = find(t>=250,1)+1:find(t>=500,1)-1; %0 to delay period
        tDTind3 = find(t>=500,1)+1:find(t>=750,1)-1; %0 to delay period
        tDTind4 = find(t>=750,1)+1:find(t>=1000,1)-1; %0 to delay period 
        tDTind5 = find(t>=1000,1)+1:find(t>=1250,1)-1; %0 to delay period 
        tDTind6 = find(t>=1250,1)+1:find(t>=1500,1)-1; %0 to delay period 
    end
    
    
    % set up dat
    for xxi = 1:2
        switch xxi
            case 1
                lDT = 450;%low delay time
                hDT = 550;%high delay time
            case 2
                lDT = 1450;
                hDT = 1550;
        end
        RTc=vRT(vDT>lDT & vDT <hDT); %RTs for assigned delay time
        sRTc=sort(RTc);
        rSRTc = sRTc(sRTc>0); %get rid of early errors
        bnRTc=round(length(rSRTc)*1/3);
        RTfc=rSRTc(bnRTc); %cutoff for fastest 1/3
        RTsc=rSRTc(bnRTc*2); %cutoff for slowest 1/3
        cT = [1:length(vDT)]'; %#ok<NBRAK>
        cO = round(length(vDT)/3);
        aa = cO;
        bb = length(vDT) - cO;
        ifastc = (vDT>lDT & vDT<hDT & vRT>0 & vRT<RTfc);
        islowc = (vDT>lDT & vDT<hDT & vRT>RTsc);
        iearlyc = (cT<=aa & vDT>lDT & vDT<hDT);
        ilatec = (cT>=bb & vDT>lDT & vDT<hDT);
        iallc = (vDT>lDT & vDT<hDT);
        
        switch xxi
            case 1
                ifast500 = ifastc;
                islow500 = islowc;
                iearly500 = iearlyc;
                ilate500 = ilatec;
                iall500 = iallc;
            case 2
                ifast1500 = ifastc;
                islow1500 = islowc;
                iearly1500 = iearlyc;
                ilate1500 = ilatec;
                iall1500 = iallc;
        end
    end
    
    disp('Performing single trial connectivity analysis...')
    
    
    %% Create adj_matrix using stPLV
    if splitTrials>1;%if size(cDat,1)>130 % split W_st into smaller files by trial number
        nSubsetTrials = round(size(cDat,3)/splitTrials);
        fprintf('splitting data into ~%d trial groups\n',nSubsetTrials);
        totalTrials = size(cDat,3);
        cDat_split={};
        for tSubset=1:splitTrials
            if tSubset<splitTrials
                cDat_split(:,:,tSubset) = {cDat(:,:,(tSubset-1)*nSubsetTrials+1:nSubsetTrials*tSubset)};
            else
                cDat_split(:,:,tSubset) = {cDat(:,:,(tSubset-1)*nSubsetTrials+1:end)};
            end
        end
    else
        cDat_split(:,:,1) = {cDat};
        nSubsetTrials = size(cDat,3);
        totalTrials = size(cDat,3);
    end
    for tSubset = 1:length(cDat_split)
        if tSubset<length(cDat_split)
            fprintf('calculating adjacency matrix for trials %d - %d of %d\n',(tSubset-1)*nSubsetTrials+1,nSubsetTrials*tSubset,totalTrials);
        else
            fprintf('calculating adjacency matrix for trials %d - %d of %d\n',(tSubset-1)*nSubsetTrials+1,totalTrials,totalTrials);
        end
        cDat = cDat_split{tSubset};
        W_st = struct;
        if twin == 0.5
            eTs = {'tBind1' 'tBind' 'tDTind1' 'tDTind2' 'tDTind3'};
        elseif twin == 0.25
            eTs = {'tBind3' 'tBind2' 'tBind1' 'tBind' 'tDTind1' 'tDTind2' 'tDTind3' 'tDTind4' 'tDTind5' 'tDTind6'};
        end
        for freq = 1:length(filtSpec)     % bandpass
            disp([filtSpec(freq).name ' frequency'])
            fSpec.range = filtSpec(freq).range;
            fSpec.order = filtSpec(freq).order;
            cstplv = stPLV(cDat,512,fSpec,twin);
            for xi = 1:length(eTs)
                eTi = eval(eTs{xi});
                disp(eTs{xi})
                for ai=1:size(cDat,3)
                    cW=threshold_proportional(squeeze(mean(cstplv(ai).plv(eTi,:,:), 'omitnan')),1);
                    A = cW;
                    [n,~] = size(A);
                    B = A'+A;
                    B(1:n+1:end)=diag(A);
                    cW = B;
                    if twin == 0.5
                        switch eTs{xi}
                            case 'tBind1'
                                Wtwin = [-1000 -500];
                                W_st(ai).precueone(freq).fRng = filtSpec(freq).range;
                                W_st(ai).precueone(freq).adj = cW;
                                W_st(ai).precueone(freq).twin = Wtwin;
                                W_st(ai).precueone(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).precueone(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                            case 'tBind'
                                Wtwin = [-500 0];
                                W_st(ai).precue(freq).fRng = filtSpec(freq).range;
                                W_st(ai).precue(freq).adj = cW;
                                W_st(ai).precue(freq).twin = Wtwin;
                                W_st(ai).precue(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).precue(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                            case 'tDTind1'
                                Wtwin = [0 500];
                                W_st(ai).DTone(freq).fRng = filtSpec(freq).range;
                                W_st(ai).DTone(freq).adj = cW;
                                W_st(ai).DTone(freq).twin = Wtwin;
                                W_st(ai).DTone(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).DTone(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                            case 'tDTind2'
                                Wtwin = [500 1000];
                                W_st(ai).DTtwo(freq).fRng = filtSpec(freq).range;
                                W_st(ai).DTtwo(freq).adj = cW;
                                W_st(ai).DTtwo(freq).twin = Wtwin;
                                W_st(ai).DTtwo(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).DTtwo(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                            case 'tDTind3'
                                Wtwin = [1000 1500];
                                W_st(ai).DTthree(freq).fRng = filtSpec(freq).range;
                                W_st(ai).DTthree(freq).adj = cW;
                                W_st(ai).DTthree(freq).twin = Wtwin;
                                W_st(ai).DTthree(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).DTthree(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                        end
                    elseif twin ==0.25
                        switch eTs{xi}
                            case 'tBind3'
                                Wtwin = [-1000 -750];
                                W_st(ai).precuethree(freq).fRng = filtSpec(freq).range;
                                W_st(ai).precuethree(freq).adj = cW;
                                W_st(ai).precuethree(freq).twin = Wtwin;
                                W_st(ai).precuethree(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).precuethree(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                            case 'tBind2'
                                Wtwin = [-750 -500];
                                W_st(ai).precuetwo(freq).fRng = filtSpec(freq).range;
                                W_st(ai).precuetwo(freq).adj = cW;
                                W_st(ai).precuetwo(freq).twin = Wtwin;
                                W_st(ai).precuetwo(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).precuetwo(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                            case 'tBind1'
                                Wtwin = [-500 -250];
                                W_st(ai).precueone(freq).fRng = filtSpec(freq).range;
                                W_st(ai).precueone(freq).adj = cW;
                                W_st(ai).precueone(freq).twin = Wtwin;
                                W_st(ai).precueone(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).precueone(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                            case 'tBind'
                                Wtwin = [-250 0];
                                W_st(ai).precue(freq).fRng = filtSpec(freq).range;
                                W_st(ai).precue(freq).adj = cW;
                                W_st(ai).precue(freq).twin = Wtwin;
                                W_st(ai).precue(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).precue(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                            case 'tDTind1'
                                Wtwin = [0 250];
                                W_st(ai).DTone(freq).fRng = filtSpec(freq).range;
                                W_st(ai).DTone(freq).adj = cW;
                                W_st(ai).DTone(freq).twin = Wtwin;
                                W_st(ai).DTone(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).DTone(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                            case 'tDTind2'
                                Wtwin = [250 500];
                                W_st(ai).DTtwo(freq).fRng = filtSpec(freq).range;
                                W_st(ai).DTtwo(freq).adj = cW;
                                W_st(ai).DTtwo(freq).twin = Wtwin;
                                W_st(ai).DTtwo(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).DTtwo(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                            case 'tDTind3'
                                Wtwin = [500 750];
                                W_st(ai).DTthree(freq).fRng = filtSpec(freq).range;
                                W_st(ai).DTthree(freq).adj = cW;
                                W_st(ai).DTthree(freq).twin = Wtwin;
                                W_st(ai).DTthree(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).DTthree(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                            case 'tDTind4'
                                Wtwin = [750 1000];
                                W_st(ai).DTfour(freq).fRng = filtSpec(freq).range;
                                W_st(ai).DTfour(freq).adj = cW;
                                W_st(ai).DTfour(freq).twin = Wtwin;
                                W_st(ai).DTfour(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).DTfour(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                            case 'tDTind5'
                                Wtwin = [1000 1250];
                                W_st(ai).DTfive(freq).fRng = filtSpec(freq).range;
                                W_st(ai).DTfive(freq).adj = cW;
                                W_st(ai).DTfive(freq).twin = Wtwin;
                                W_st(ai).DTfive(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).DTfive(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                            case 'tDTind6'
                                Wtwin = [1250 1500];
                                W_st(ai).DTsix(freq).fRng = filtSpec(freq).range;
                                W_st(ai).DTsix(freq).adj = cW;
                                W_st(ai).DTsix(freq).twin = Wtwin;
                                W_st(ai).DTsix(freq).trial = ai+((tSubset-1)*nSubsetTrials);
                                W_st(ai).DTsix(freq).RT = vRT(ai+((tSubset-1)*nSubsetTrials));
                        end
                    end
                end
            end
            clear cstplv
        end
        
        %% Save
        if ~exist(saveDir,'dir'),mkdir(saveDir);end
        if saveon
            disp('Saving...') %#ok<UNRCH>
            if length(cDat_split)>1
                save(fullfile(saveDir,[subj '_Wst' num2str(bipol) '_' num2str(tSubset)]), 'W_st', 't*', 'cDat', 'ch*', 'JAC', 'vRT', 'vDT', 'i*');
            else
                save(fullfile(saveDir,[subj '_Wst' num2str(bipol)]), 'W_st', 't*', 'cDat', 'ch*', 'JAC', 'vRT', 'vDT', 'i*');
            end
        end
    end
    disp('Done');
end