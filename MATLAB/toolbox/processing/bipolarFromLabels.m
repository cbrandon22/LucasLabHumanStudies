function bpCells = bipolarFromLabels(labels,channels)
%%% This function creates bipolar pairs based only on labels. If you
%%% also provide matched channel numbers, that will be included as a second
%%% column in the output
[~,~,mapCells]=xlsread('electrodeMap.xlsx');
channels = mapCells(1:end-1,1);
labels = mapCells(1:end-1,3);
clear mapCells

if exist('channels','var')
    if ~isequal(size(channels),size(labels))
        error('Labels and channels dimensions must match');
    end
end
cellfun('@(s)' strsplit(s,'\d','DelimiterType','RegularExpression'),labels)
regexp(labels,'\d')
[lbl,matches] = strsplit(labels{1},'\d','DelimiterType','RegularExpression');
lbl=lbl{1};eNum=matches{1};
while
