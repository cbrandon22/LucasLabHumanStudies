function powPlot
subj = 'HUP119_i';
dirs = an_dirs;
powerDir = fullfile(dirs.scratch,subj);
files = dir(fullfile(powerDir,'*.mat'));
for j=1: length(files)
    file = files(j).name;
    tank = strsplit(file,[subj '_']); %FIX regexp***
    load(file);
    for ich = 1:size(PM,1)
        PM(ich,:,:) = smooth2(squeeze(PM(ich,:,:)),2,1);
        PM(ich,:,:) = zscore(squeeze(PM(ich,:,:))')'; % zscore
    end
    
    % plot
    for ich = 1:size(PM,1)
        figure('Name',[subj ' Electrode: ' tank],'NumberTitle','off','Units','normalized','Position',[1/8 1/4 3/4 1/2],'Color','w');
        imagesc(x,1:length(freq),squeeze(PM(ich,:,:)));
        ytick = [1:10 20:10:100 200];
        yticknm = {'1';'';'';'';'';'';'';'';'';'10';'';'';'';'';'';'';'';'';'100';''};
        yind = zeros(size(ytick));
        for ii = 1:length(ytick)
            [~,yind(ii)] = min(abs(freq-ytick(ii)));
        end
        set(gca,'Box','off','YDir','normal','TickDir','out','YLim',[1 length(freq)],'XLim',[x(1) x(end)],'YTick',yind,'YTickLabel',yticknm);
        xlabel('time (min)'); ylabel('frequency (Hz)'); title([subj ' Electrode: ' tank],'Interpreter','none');
        set(gca,'CLim',[-3 3]);
        hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
    end
end

end