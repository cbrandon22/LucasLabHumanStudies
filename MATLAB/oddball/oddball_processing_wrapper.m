% subjList = {'HUP142_i','HUP144_e','HUP145_e','HUP147_e',...
%     'HUP148_e','HUP149_e','HUP150_i','HUP151_e','HUP152_e','HUP153_i','HUP154_e',...
%     'HUP155_i','HUP156_i','HUP157_e','HUP159_e','HUP166_i'};
subjList = {'HUP166_i'};
for i=1:length(subjList)
    disp(subjList{i})
    an_processSubj2(subjList{i},0);
    oddball_processing(subjList{i});
end