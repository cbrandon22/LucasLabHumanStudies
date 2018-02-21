function [chi2stat,p,df,O,E]= doChi2Stat(id,expFreq,totCount)
%This is a wrapper for chi2gof computes chiSq stats.
%Inputs:
%id                     %vector categorizing each observation (e.g, 0s and 1)
%expFreq = [.5 .5]      %expected frequencies associated with each category     
%totCount (optional)    %total num of observations. Enter this value if
                        %length(id)>totCount
 
%Outputs:
%chi2stat
%p
%df
%Observed
%Expected

%Written by AGR 03-03-13
%if only one value in idVec return
if length(unique(id))==1
    chi2stat = nan; p = nan; df = nan; O = [length(id) 0]; E = expFreq*length(id);
    return
end


if ~exist('totCount','var') || isempty(totCount)
    totCount = length(id);
end
if totCount<length(id)
    id(find(id==0,length(id)-totCount)) = [];
end
% it doesn't like 0's so we increment all values by 1
if sum(id==0)>0 
    idVec = id + 1;
else
    idVec = id;
end
uniqueVals = unique(idVec);

[~,p,stats]=...
    chi2gof(idVec,'ctrs',uniqueVals,'expected',[expFreq*totCount],'emin',0);
chi2stat = stats.chi2stat;
df = stats.df;
O = stats.O;
E = stats.E;
