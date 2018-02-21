function [h] = swagFig(pos)
if ~exist('pos','var') || isempty(pos)
    pos = [0.7504    0.6850    0.2188    0.2662];
end
h = figure('units','normalized','paperpositionmode','auto','position',pos);
