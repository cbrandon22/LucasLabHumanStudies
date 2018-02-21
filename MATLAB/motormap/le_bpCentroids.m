%Save .csv file of bipolar centroids
subj = 'HUP089';

dirs = le_dirs;
dataDir = dirs.data;

[anatStruct_bipol] = le_centroids2Anatstruct(subj);
centroids = [];
for i = 1:length(anatStruct_bipol)
    centroids(i,:) = [anatStruct_bipol(i).X, anatStruct_bipol(i).Y, anatStruct_bipol(i).Z];
end

csvwrite(fullfile(dataDir,'eeg',subj,'tal','bpCentroids.csv'),centroids);