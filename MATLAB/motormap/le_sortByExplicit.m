function [sortIdx,exp_feat] = sortByExplicit(tMat_lf_cat,...
    tStruct_lf,config_mm,config_pow)
fBins = [config_pow.freqBins];
fBins = fBins(tStruct_lf(1).fInd); % filter by freq range used for this tStruct


switch config_mm.explicit_feature_to_sort_by
    case {'theta','beta'}
        if strcmp(config_mm.explicit_feature_to_sort_by,'theta')
            fRange = config_mm.featfRange.theta;
        elseif strcmp(config_mm.explicit_feature_to_sort_by,'beta')
            fRange = config_mm.featfRange.theta;
        end
        [~,fInd_start] = min(abs(fRange(1) - fBins));
        [~,fInd_end] = min(abs(fRange(2) - fBins));
        fInd = fInd_start:fInd_end;

        % collapse based on that feature
        exp_feat = nanmean(tMat_lf_cat(:,fInd),2);
        [~,sortIdx] = sort(exp_feat);

case {'diff(theta,beta)'}
        %collapse theta pow
        fRange = config_mm.featfRange.theta;
        [~,fInd_start] = min(abs(fRange(1) - fBins));
        [~,fInd_end] = min(abs(fRange(2) - fBins));
        fInd = fInd_start:fInd_end;
        theta_pow = nanmean(tMat_lf_cat(:,fInd),2);
        
        %collapse beta pow
        fRange = config_mm.featfRange.beta;
        [~,fInd_start] = min(abs(fRange(1) - fBins));
        [~,fInd_end] = min(abs(fRange(2) - fBins));
        fInd = fInd_start:fInd_end;
        beta_pow = nanmean(tMat_lf_cat(:,fInd),2);
        
        exp_feat = [theta_pow beta_pow];
        [~,sortIdx] = sort(exp_feat(:,1) - exp_feat(:,2));

end