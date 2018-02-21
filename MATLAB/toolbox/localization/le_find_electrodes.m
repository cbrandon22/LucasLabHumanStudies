%%
%This script can be used to find centroids for electrode locations that
%have been generated using brainmapper and ITK-SNAP.

dirs = le_dirs;
dataDir = dirs.data;
addpath('/Volumes/Lucas_ECoG_Data');
% USER INPUTS
subj = 'HUP089';

%segmentation label intensity
label = 2;
seg = load_untouch_nii('/Volumes/Lucas_ECoG_Data/data/eeg/HUP089/images/HUP089_unburied_electrode_seg.nii');
%%
untouched_img = seg.img;
electrode_map = single(untouched_img>=label);
CC = bwconncomp(electrode_map);
S = regionprops(CC,'Centroid');
electrode_vis = electrode_map + single(untouched_img>0);
centroids = cat(1,S.Centroid);

% swaps x and y axis values to match up with original image
centroids(:,[1,2]) = centroids(:,[2,1]);
%%
csvwrite(fullfile(dataDir,'eeg',subj,'tal','centroids.csv'),centroids);
