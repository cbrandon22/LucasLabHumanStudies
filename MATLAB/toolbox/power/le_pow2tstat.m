
function[t,p,pow1,pow2,fInd,tInd] = le_pow2tstat(powMat,retain1,retain2,fRange,tRange,collapseFreqFlag,config);
%outputs
%t = vector of t's
%p ... vector of p's
%pow1... power vector for condition 1 
%pow2 ... power vector for condition 2
%fInd .... freqs of interest
%tInd ... time points of interest
 %this computes a t-stat and collapses across fRange and tRange
%keyboard 

%get freq range
[~,fInd_start] = min(abs(fRange(1) - config.freQ));
[~,fInd_end] = min(abs(fRange(2) - config.freQ));
fInd = fInd_start:fInd_end;

% get time range
[~,tInd_start] = min(abs(tRange(1) - nanmean(config.timeBins,2)));
[~,tInd_end] = min(abs(tRange(2) - nanmean(config.timeBins,2)));
tInd = tInd_start:tInd_end;


%compute power vectors and do ttests (separate for collapseFreqFlag =1 and 0)
if collapseFreqFlag % TRUE (this collapses across the freq range)
    
    pow1 =  shiftdim(nanmean(powMat(fInd,tInd,retain1),1),1)';
    if sum(retain2)>0
        pow2 =  shiftdim(nanmean(powMat(fInd,tInd,retain2),1),1)';
    else
        pow2 = [];
    end

    % do unpaired t-test
    if isempty(retain2) || sum(retain2) == 0 % one sample
        [~,p,~,stats] = ttest(pow1,0,'dim',1);
        t = stats.tstat;
    else % do a two sample
        [~,p,~,stats] = ttest2(pow1,pow2,'dim',1);
        t = stats.tstat;
    end

else % FALSE (this performs a time frequency analysis)
    
    pow1 =  powMat(fInd,tInd,retain1);
    if sum(retain2)>0
        pow2 =  powMat(fInd,tInd,retain2);
    else
        pow2 = [];
    end

    % do unpaired t-test
    if isempty(retain2) || sum(retain2) == 0 % one sample
        [~,p,~,stats] = ttest(pow1,0,'dim',3);
        t = stats.tstat;
    else % do a two sample
        [~,p,~,stats] = ttest2(pow1,pow2,'dim',3);
        t = stats.tstat;
    end

    
    
end