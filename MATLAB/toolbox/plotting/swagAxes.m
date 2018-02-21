function[]=swagAxes(a,fsize,xlbl,ylbl,tit,loose)
%this function takes a figure and swags it out / makes it look cool
%Inputs
%a.........axes handle
%optional
%xlbl.........xlabel ('str')
%ylbl.........ylabel ('str')
%title.........title ('str')
%loose.........if given a value, it provides a little space around the
               %edges of the plotted data 

%a=axes handle
axis('tight')
set(a,'fontsize',fsize,'box','off')
if exist('xlbl','var') && ~isempty(xlbl)
xlabel(xlbl,'fontsize',fsize)
end
if exist('ylbl','var') && ~isempty(ylbl)
ylabel(ylbl,'fontsize',fsize)
end
if exist('tit','var') && ~isempty(tit)
title(tit,'fontsize',(fsize*2))
end

if exist('loose','var') && ~isempty(loose)
%loosen the axes a little bit
xl = get(gca,'xlim');
xRange = xl(2) - xl(1);
set(gca,'xlim',[xl(1)-(.05*xRange) xl(2)+(.05*xRange)])

yl = get(gca,'ylim');
yRange = yl(2) - yl(1);
set(gca,'ylim',[yl(1)-(.05*yRange) yl(2)+(.05*yRange)])
end