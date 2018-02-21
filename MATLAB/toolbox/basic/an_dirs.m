function [dirs] = an_dirs()
%function maintains a list of directories for ecog analyses
root = '/Volumes/Lucas_ECoG_Data/Consciousness';

dirs.data = fullfile(root,'data');
dirs.scratch = fullfile(root,'scratch');
dirs.lockDir = fullfile(root,'lockDir');
dirs.reportDir = fullfile(root,'ecogReports');
