function ram_coh_interarray
% function ram_coh_interarray
%   Plot matrix of coherence spectra between pairs of recording sites on
%   two different 8-electrode arrays.
% 
%   DR 06/2015

% parameters
ddir = 'H:\ditka\'; % directory
dtank = 'ditka_rec_DT1_033015'; % data tank
dblock = 1; % data block
a1 = 1:8; % channels: 1st array
a2 = 9:16; % channels: 2nd array
freq = logspace(log10(0.3),log10(50),25); % frequencies (Hz)
dns = 10; % downsample factor

% data
cd([ddir dtank '\Block-' num2str(dblock)]);
fs = 24414.06; % sampling rate (Hz)
dfs = round(fs/dns);
[bA,aA] = butter(2,500/(fs/2)); % antialiasing
figure('Name',dtank,'NumberTitle','off','Units','normalized','Position',[1/4 1/8 1/2 3/4],'Color','w');
x = linspace(-0.4,0.4,length(freq));
for ii = a1
    fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(ii) '.sev'],'r');
    fseek(fid,40,'bof');
    dat1 = fread(fid,inf,'single');
    fclose(fid);
    dat1 = filtfilt(bA,aA,dat1);
    dat1 = downsample(dat1',dns)';
    for jj = a2
        fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(jj) '.sev'],'r');
        fseek(fid,40,'bof');
        dat2 = fread(fid,inf,'single');
        fclose(fid);
        dat2 = filtfilt(bA,aA,dat2);
        dat2 = downsample(dat2',dns)';

%         Cxy = mscohere(dat1,dat2,10*dfs,5*dfs,freq,dfs);
%         plot(x+jj-a2(1),Cxy+ii-a1(1),'k'); hold on;
        [Cxy,f] = mscohere(dat1,dat2,10*dfs,5*dfs,10*dfs,dfs);
        Cxyf = zeros(size(freq));
        for kk = 1:length(freq)
            [~,ind] = min(abs(freq(kk)-f));
            Cxyf(kk) = Cxy(ind);
        end
        plot(x+jj-a2(1),Cxyf+ii-a1(1),'k'); hold on;
        
        pause(1);
    end
end
axis square; set(gca,'Box','off','XLim',[-1 length(a2)],'XTick',0:length(a2)-1,'XTickLabel',a2,'YLim',[0 length(a1)],'YTick',0.5:length(a1),'YTickLabel',a1);
xlabel('channel'); ylabel('channel');