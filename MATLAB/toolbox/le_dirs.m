function [dirs] = le_dirs(task)
%function maintains a list of directories for ecog analyses
root = fullfile('/Volumes/HumanStudies/HumanStudies',task);
if strcmp(task,'motormap')
    root = fullfile('/Volumes/HumanStudies 1/HumanStudies',task);
end

dirs.data = fullfile(root);
dirs.events = fullfile(root,'events');
dirs.scratch = fullfile(root,'scratch');
dirs.lockDir = fullfile(root,'lockDir');
dirs.reportDir = fullfile(root,'ecogReports');