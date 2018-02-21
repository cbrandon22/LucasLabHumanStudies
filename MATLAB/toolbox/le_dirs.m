function [dirs] = le_dirs(task)
%function maintains a list of directories for ecog analyses
root = fullfile('/Volumes/HumanStudies/HumanStudies',task);

dirs.data = fullfile(root);
dirs.events = fullfile(root,'events');
dirs.scratch = fullfile(root,'scratch');
dirs.lockDir = fullfile(root,'lockDir');
dirs.reportDir = fullfile(root,'ecogReports');