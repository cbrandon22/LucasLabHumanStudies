%%
% This script will take one centroid coordinate, create a spherical ROI (of
% 'radius' voxels in radius) around that one centroid and output a
% segmentation file that can be loaded into ITK-SNAP.

%% Load external libraries
addpath('~/Documents/MATLAB/NIfTI_20140122');

%% Parameters
dirs = le_dirs;
dataDir = dirs.data;

% USER INPUTS
subj = 'HUP111';

radius = 5; %voxels

%% Load centroids as well as seg
% Load segmentation label intensity
label = 2;
seg = load_untouch_nii(sprintf('/Volumes/Lucas_ECoG_Data/data/eeg/%s/images/ReconstructionSegImg/%s_2nd_unburied_electrode_seg.nii.gz',subj,subj));
seg = seg.img;

% Brain Segmentation
brain_mask = zeros(size(seg));
brain_mask(seg>=1) = 1;

% Load centroids
centroids = csvread(fullfile(dataDir,'eeg',subj,'tal','centroids.csv'));
% TODO: Figure out how to cleanly load labels.csv into table/array
% labels = csvread(fullfile(dataDir,'eeg',subj,'tal','labels.csv'));

%% Pick one centroid 
centroid_id = 10;
ele_centroid = centroids(10,:);

% ROI mask
ele_roi = zeros(size(seg));
ele_roi(brain_mask == 1) = 1;
ele_roi(ele_centroid(1)-radius:ele_centroid(1)+radius,ele_centroid(2)-radius:ele_centroid(2)+radius,ele_centroid(3)-radius:ele_centroid(3)+radius) = 2;
ele_roi(brain_mask == 0) = 0;


%% Save new electrode ROI
seg = load_untouch_nii(sprintf('/Volumes/Lucas_ECoG_Data/data/eeg/%s/images/ReconstructionSegImg/%s_2nd_unburied_electrode_seg.nii.gz',subj,subj));
seg.img = ele_roi;
save_untouch_nii(seg,sprintf('/Volumes/Lucas_ECoG_Data/data/eeg/%s/tal/test_electrode_roi_LAST4.nii.gz',subj));