% This script builds group powDiff struct and collapses anatomical regions
% for plot.

% load powDiff files for each subject
subjList = {'HUP084','HUP087','HUP091','HUP092','HUP111'};
[powDiffBatch] = vb_batchLoad(subjList,'_powDiff.mat',{'powDiff'});

% build group structure
powDiffStruct = struct();
for i=1:length(powDiffBatch)
    if i==1
        for j=length(powDiffStruct):length(powDiffBatch{i}.powDiff)
            powDiffStruct(j).subj = {powDiffBatch{i}.name};
            powDiffStruct(j).pow = {powDiffBatch{i}.powDiff(j).pow};
            powDiffStruct(j).eLbl = [powDiffBatch{i}.powDiff(j).eLbl];
            powDiffStruct(j).anatAbbr = [powDiffBatch{i}.powDiff(j).anatAbbr];     
        end
    else
        startInd = length(powDiffStruct);
        for j=(1+startInd):(length(powDiffBatch{i}.powDiff)+startInd)
            powDiffStruct(j).subj = {powDiffBatch{i}.name};
            powDiffStruct(j).pow = {powDiffBatch{i}.powDiff(j-startInd).pow};
            powDiffStruct(j).eLbl = [powDiffBatch{i}.powDiff(j-startInd).eLbl];
            powDiffStruct(j).anatAbbr = [powDiffBatch{i}.powDiff(j-startInd).anatAbbr];     
        end
    end
end

% Sort powDiffStruct by anatAbbr
[tmp,ind]=sort([powDiffStruct.anatAbbr]);
powDiffStruct = powDiffStruct(ind);

% Create list of each anatomical region included
anatList = {};
for i=1: length(powDiffStruct)
   if ~ismember(powDiffStruct(i).anatAbbr,anatList)
       anatList(length(anatList)+1) = powDiffStruct(i).anatAbbr;
   end
end

% Collapse to one powDiff vector per anatomical region
anatPowStruct = struct();
j=1;
for i=1: length(anatList)
    % build tempStruct for each anatomical location
    tempStruct = struct();
    x=1;
   while (j<=length(powDiffStruct) && strcmp(anatList{i}, powDiffStruct(j).anatAbbr))
       tempStruct(x).anatAbbr = anatList{i};
       tempStruct(x).pow = powDiffStruct(j).pow;
       j = j+1;
       x = x+1;
   end
   % collapse tempstruct to pow vector across frequencies and add to
   % anatPowStruct
   anatPowStruct(i).anatAbbr = anatList{i};
   anatPowStruct(i).n = length(tempStruct);
   if length(tempStruct)==1
       anatPowStruct(i).pow = squeeze(nanmean(tempStruct.pow{1},2));
   else
       powCat = tempStruct(1).pow{1};
       for y=2:length(tempStruct)
           powCat = cat(2,powCat,tempStruct(y).pow{1});
       end
       anatPowStruct(i).pow = squeeze(nanmean(powCat,2));
   end
end

powVectors = anatPowStruct(1).pow;
for i=2:length(anatPowStruct)
   powVectors = cat(2,powVectors,anatPowStruct(i).pow);
end
% plot
figure;
pcolor([1:length(anatPowStruct(1).pow)],[1:length(anatPowStruct)], (real((powVectors)))');
colorbar;
shading flat
caxis ([-0.5 0.5])
yt = [1:1:length(anatList)];
set(gca,'ytick',yt,'yticklabel',anatList(yt),'yticklabelrotation',0)

[h,p] = ttest(powVectors)
