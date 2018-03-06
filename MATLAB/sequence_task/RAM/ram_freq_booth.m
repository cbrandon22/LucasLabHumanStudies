function ram_freq_booth
% function ram_freq_booth
%   Compute power spectrum over entire session collectd in booth.
% 
%   DR 03/2015

% parameters
ddir = 'G:\'; % directory
subj = 'ditka'; % subject
tank = 'ditka_stim_DT1_032715'; % tank 
ch = 4; % channel

% main
fs = 24414.06;
try
    load([ddir subj '\processed\' tank '.mat'],'-mat');
    eval(['P = Pch' num2str(ch) ';']);
catch
    cd([ddir subj '\' tank]);
    blocks = dir('Block-*');
    bnm = zeros(length(blocks),1);
    for ii = 1:length(blocks)
        nm = blocks(ii).name;
        idash = strfind(nm,'-');
        bnm(ii) = str2double(nm(idash+1:end));
    end
    [~,is] = sort(bnm);
    freq = logspace(log10(2),log10(100),100)'; % frequencies (Hz)
    P = []; indb = [];
    for iblock = 1:length(blocks)
        indb = [indb size(P,2)+1];
        cd([ddir subj '\' tank '\' blocks(is(iblock)).name]);
        fid = fopen([tank '_' blocks(is(iblock)).name '_xRAW_ch25.sev'],'r');
        fseek(fid,40,'bof');
        dat = fread(fid,inf,'single');
        fclose(fid);
        itrain = find(dat(1:end-1)<0.5 & dat(2:end)>0.5); % onset of each pulse train (samples)
        itrain = itrain + round(0.5/1000*fs);
        itrain = [itrain; itrain + round(.5*fs)]; % onsets and offsets of each pulse train period
        
        % data
        fid = fopen([tank '_' blocks(is(iblock)).name '_xRAW_ch' num2str(ch) '.sev'],'r');
        fseek(fid,40,'bof');
        dat = fread(fid,inf,'single');
        fclose(fid);
        N = length(dat); % samples
        S = fix(N/fs); % seconds
        
        % spectra
        PS = zeros(length(freq),S);
        ss = zeros(1,S);
        fprintf('computing spectra...');
        for ii = 1:S % compute in 1-sec windows in 1-sec steps
            fprintf('%02d%%',fix(ii/S*100));
            ind = round((ii-1)*fs+1:ii*fs);
            PS(:,ii) = periodogram(dat(ind),hamming(length(ind)),freq,fs);
            if any(intersect(ind,itrain)), ss(ii) = 1; end % keep track of which seconds included stim artifact
            if ii~=S, fprintf('\b\b\b'); end
        end
        fprintf('\n');
        
        % compile
        PS(:,ss==1) = [];
        P = [P PS];
    end
    
    % save
    eval(['Pch' num2str(ch) ' = P;']);
    try eval(['save([ddir subj ''\processed\'' tank],''Pch' num2str(ch) ''',''freq'',''indb'',''-append'');']);
    catch eval(['save([ddir subj ''\processed\'' tank],''Pch' num2str(ch) ''',''freq'',''indb'');']); end
end

% plot 1: spectrogram
Pt = smooth2(P,45,2); % sec, Hz 
Pt = zscore(Pt')';
figure('Name',tank,'NumberTitle','off','Units','normalized','Position',[1/8 1/4 3/4 1/2],'Color','w');
imagesc(1:size(Pt,2)/60,1:size(Pt,1),Pt); hold on;
ytick = [2:10 20:10:100];
yticknm = {'';'';'';'';'';'';'';'';'10';'';'';'';'';'';'';'';'';'100'};
yind = zeros(size(ytick));
for ii = 1:length(ytick)
    [~,yind(ii)] = min(abs(freq-ytick(ii)));
end
set(gca,'Box','off','YDir','normal','TickDir','out','YLim',[1 length(freq)],'YTick',yind,'YTickLabel',yticknm);
xlabel('time (min)'); ylabel('frequency (Hz)'); title(tank,'Interpreter','none');
hc = colorbar; set(get(hc,'YLabel'),'String','power');
for ii = 2:length(indb)
    line(indb(ii)/60*ones(1,2),get(gca,'YLim'),'Color','w','LineWidth',2);
end

% plot 2: compare beginning to end
N = 600; % sec to compare
P1 = P(:,1:N)';
P2 = P(:,end-N+1:end)';
figure('Name',tank,'NumberTitle','off','Units','normalized','Position',[1/4 1/8 1/2 1/2],'Color','w');
m1 = mean(P1); v1 = 1.96*std(P1)/sqrt(N);
m2 = mean(P2); v2 = 1.96*std(P2)/sqrt(N);
patch([freq', fliplr(freq')],[m1-v1, fliplr(m1+v1)],'k','FaceColor',[0.8 0.8 1],'EdgeColor','none'); hold on;
patch([freq', fliplr(freq')],[m2-v2, fliplr(m2+v2)],'k','FaceColor',[1 0.8 0.8],'EdgeColor','none'); hold on;
h(1) = plot(freq,m1,'Color',[0 0 1],'LineWidth',2);
h(2) = plot(freq,m2,'Color',[1 0 0],'LineWidth',2);
axis tight; axis square; set(gca,'Box','off','XScale','log','YScale','log');
xlabel('Hz'); ylabel('\muV^2/Hz');
legend(h,['first ' num2str(N) ' sec'],['last ' num2str(N) ' sec']);
legend('boxoff');