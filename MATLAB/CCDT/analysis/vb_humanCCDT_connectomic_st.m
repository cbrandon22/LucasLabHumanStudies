% CCDT Human Connectomic Analysis
for iji = 6 % number of subjects
    subj = allSubj{iji}; %subject name string vector
    disp(subj)
    saveon=0;
    bipol=0; % 1 if want to use bipolar montage, 0 if monopolar
    lpc=200;
    fs=512;
    twin=0.5; % milliseconds for single-trial PLV sliding window
    
    %%
    disp('load events...')
    load(['D:\CCDT\events\' subj '_events.mat']);
    
    fid = fopen(['D:\CCDT\eeg\' subj '\docs\jacksheet.txt']);
    C = textscan(fid,'%s');
    fclose(fid);
    JAC{:,1} = C{1}(1:2:end-1);
    JAC{:,2} = C{1}(2:2:end);
    
    durationMS=3000;
    offsetMS=-1000;
    bufferMS=1000;
    resampleFreq=512;
    chrec = cellfun(@str2double, JAC{1});
    chlbl = JAC{2};
    
    
    bT = .300; %(s) buffer time after end of DT
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
        disp(['ch ' num2str(chrec(jj))]);
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
    
    
    tBind = find(t>=-500,1)+1:find(t>=0,1)-1; %-500 to 0ms baseline
    tDTind1 = find(t>=0,1)+1:find(t>=500,1)-1; %0 to delay period
    tDTind2 = find(t>=500,1)+1:find(t>=1000,1)-1; %0 to delay period
    tDTind3 = find(t>=1000,1)+1:find(t>=1500,1)-1; %0 to delay period    
    
    
    % set up dat
    for xxi = 1:2
        switch xxi
            case 1
                lDT = 450;
                hDT = 550;
            case 2
                lDT = 1450;
                hDT = 1550;
        end
        RTc=vRT(vDT>lDT & vDT <hDT);
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
    
    fT.range=[3 12];
    fT.order=250;
    fG.range=[70 100];
    fG.order=50;
    
    disp('Performing single trial connectivity analysis...')
    
    
    %% Create adj_matrix using stPLV
    W_st = struct;
    eTs = {'tBind' 'tDTind1' 'tDTind2' 'tDTind3'};
    for freq = 1:2     % bandpass
        switch freq
            case 1 %theta
                [cstplv] = stPLV(cDat,512,fT,twin);
            case 2 %high gamma
                [cstplv] = stPLV(cDat,512,fG,twin);
        end
        for xi = 1:length(eTs)
            eTi = eval(eTs{xi});
            for ai=1:size(cDat,3)
                cW=threshold_proportional(squeeze(mean(cstplv(ai).plv(eTi,:,:), 'omitnan')),1);
                A = cW;
                [n,~] = size(A);
                B = A'+A;
                B(1:n+1:end)=diag(A);
                cW = B;
                
                disp(eTs{xi})
                switch eTs{xi}
                    case 'tBind'
                        Wtwin = [-500 0];
                        switch freq
                            case 1 %theta
                                W_st(ai).precue.theta.adj = cW;
                                W_st(ai).precue.theta.twin = Wtwin;
                                W_st(ai).precue.theta.trial = ai;
                                W_st(ai).precue.theta.RT = vRT(ai);
                            case 2 %high gamma
                                W_st(ai).precue.gamma.adj = cW;
                                W_st(ai).precue.gamma.twin = Wtwin;
                                W_st(ai).precue.gamma.trial = ai;
                                W_st(ai).precue.gamma.RT = vRT(ai);
                        end
                    case 'tDTind1'
                        Wtwin = [0 500];
                        switch freq
                            case 1 %theta
                                W_st(ai).DTone.theta.adj = cW;
                                W_st(ai).DTone.theta.twin = Wtwin;
                                W_st(ai).DTone.theta.trial = ai;
                                W_st(ai).DTone.theta.RT = vRT(ai);
                            case 2 %high gamma
                                W_st(ai).DTone.gamma.adj = cW;
                                W_st(ai).DTone.gamma.twin = Wtwin;
                                W_st(ai).DTone.gamma.trial = ai;
                                W_st(ai).DTone.gamma.RT = vRT(ai);
                        end
                    case 'tDTind2'
                        Wtwin = [500 1000];
                        switch freq
                            case 1 %theta
                                W_st(ai).DTtwo.theta.adj = cW;
                                W_st(ai).DTtwo.theta.twin = Wtwin;
                                W_st(ai).DTtwo.theta.trial = ai;
                                W_st(ai).DTtwo.theta.RT = vRT(ai);
                            case 2 %high gamma
                                W_st(ai).DTtwo.gamma.adj = cW;
                                W_st(ai).DTtwo.gamma.twin = Wtwin;
                                W_st(ai).DTtwo.gamma.trial = ai;
                                W_st(ai).DTtwo.gamma.RT = vRT(ai);
                        end
                    case 'tDTind3'
                        Wtwin = [1000 1500];
                        switch freq
                            case 1 %theta
                                W_st(ai).DTthree.theta.adj = cW;
                                W_st(ai).DTthree.theta.twin = Wtwin;
                                W_st(ai).DTthree.theta.trial = ai;
                                W_st(ai).DTthree.theta.RT = vRT(ai);
                            case 2 %high gamma
                                W_st(ai).DTthree.gamma.adj = cW;
                                W_st(ai).DTthree.gamma.twin = Wtwin;
                                W_st(ai).DTthree.gamma.trial = ai;
                                W_st(ai).DTthree.gamma.RT = vRT(ai);
                        end
                end
            end
        end
        clear cstplv
    end
    
    disp('Done');
    
    %% Save
    
    if saveon
        disp('Saving...') %#ok<UNRCH>
        save([subj '_Wst' num2str(bipol)], 'W_st', 't*', 'cDat', 'ch*', 'JAC', 'vRT', 'vDT', 'i*');
        keep allSubj
    end
        
end