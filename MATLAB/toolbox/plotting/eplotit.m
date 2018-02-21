function out = eplotit(dat,t,trialname,yl,z)
% dat is a trial x time matrix for each comparison of data desired.
    % e.g. out = eplotit({Pca{1} Pia{1} PiaP{1}},t,'Correct vs. Error vs. Sedation')
% t = time
% edh 2017 ehudgins@gmail.com
%
%
% color = '-k'; % line color, 
% DAT color scheme goes blue = 1, red = 2, yellow = 3
if isempty(z)
    z = 0;
end
if isempty(yl)
    yl = 'Average Power';
end
fillcol = [0.7 0.7 1]; % fill color
cl4 = fillcol; % [0.7 0.7 1]; % fill color
transp = 0.8; % transparency for fill

for n = 1:length(dat) % get number of comparisons to add to the plot
    out.M{n} = nanmean(dat{n});
    out.S{n} = nanstd(dat{n});
    out.E{n} = out.S{n}/sqrt(size(dat{n},1));
    if n < length(dat)
        [~,out.Pval{n}] = ttest2(dat{n},dat{n+1},'tail','both');
        [out.pthr{n},out.pcor{n},out.padj{n}] = fdr2(out.Pval{n});
%         [out.pID{n},out.pN{n}] = fdr(out.Pval{n},0.05);
    else
        [~,out.Pval{n}] = ttest2(dat{n},dat{1},'tail','both');
        [out.pthr{n},out.pcor{n},out.padj{n}] = fdr2(out.Pval{n});
%         [out.pID{n},out.pN{n}] = fdr(out.Pval{n},0.05);
    end
end

%  fill([out.t1 out.t1(end:-1:1)], [out.y1 out.y2(end:-1:1)], 'w','EdgeColor', cl4, 'EdgeAlpha',transp, 'FaceColor', cl4, 'FaceAlpha',transp)
%  plot(out.t1, out.mr, color, 'linewidth', 1.5)
 
 figure('name',[trialname ' mean SEM'])
 for n = 1:length(out.M)
     if z == 1
         M = zscore(out.M{n});
%          E = std(zscore(dat{n}))/sqrt(size(dat{n},1));
         E = zscore(out.E{n});
     else
         M = out.M{n};
         E = out.E{n};
     end
      y1 = M + E;
      y2 = M - E;
     fill([t t(end:-1:1)], [y1 y2(end:-1:1)], 'w','EdgeColor', cl4, 'EdgeAlpha',transp, 'FaceColor', cl4, 'FaceAlpha',transp)
     hold on
     h=plot(t, M, 'linewidth', 2);
     plot(t(find(out.Pval{n} < 0.0001)),M(find(out.Pval{n} < 0.0001)),'+y');
     out.h{n}=h;
 end
    xlabel('time (ms)'); ylabel(yl); title(trialname);
    line([0 0],get(gca,'YLim'),'Color','k');
    

% Early = Pa(1:round(size(Pa,1)/2),:); mEarly = mean(Early); sEarly = std(Early); eEarly = sEarly/sqrt(size(Early,1));
% Late = Pa(round(size(Pa,1)/2)+1:size(Pa,1),:); mLate = mean(Late); sLate = std(Late); eLate = sLate/sqrt(size(Late,1));
% 
% figure('name',[trialname 'mean STdev'])
% errorbar(t,mEarly,eEarly)
% hold on
% plot(t,mEarly)
% errorbar(t,mLate,eLate)
% plot(t,mLate)
% errorbar(t,Mbl,Sbl)
% plot(t,Mbl)
% xlabel('time (ms)'); ylabel('Average Power'); title([trialname ' Hippocampus']);

