function ram_freq_cage_ditka
% function ram_freq_cage
%   Compute power spectrum from data collected wirelessly in the cage.
% 
%   DR 03/2015

% parameters
ddir = '/Users/tnl/Desktop/'; % directory
% load('/Users/tnl/Desktop/Sleep/subj.mat')
% subj.ditka = {'ditka_cage_040615' 'ditka_cage_050515' 'ditka_cage_051115' 'ditka_cage_051315' 'ditka_cage_051815' 'ditka_cage_052615' 'ditka_cage_060215'};
% subj.fidel = {'fidel_cage_082715' 'fidel_cage_093015' 'fidel_cage_100515' 'fidel_cage_100815' 'fidel_cage_101415' 'fidel_cage_111715'};; % subject

dday = [40615 50515 51115 51315 51815 52615 60215];
srates =[1000 1000 10000 1000 10000 10000 10000];
chns = [1 1 1 1 1 1 1];
subj = 'ditka';

for j=1:length(dday);

    tank = ['ditka_cage_0' num2str(dday(j))]; % tank 
    fs = srates(j); % sampling rate (Hz)
    nch = chns(j); % number of channels

% data
        freq = logspace(log10(0.5),log10(200),100)'; % frequencies (Hz)
        winsz = 5; % window size for spectral estimates (s) - determines frequency resolution
        
        try load([ddir subj filesep 'processed' filesep tank '.mat'],'-mat');
        catch
            cd([ddir subj filesep tank]);
            disp('loading data...');
            Nfl = length(dir('*_*.txt'));
            DAT = [];
                for ii = 0:Nfl-1
                    fnm = [tank(end-5:end-2) '_' num2str(ii) '.txt'];
                    fid = fopen(fnm,'r');
                    DAT = [DAT; fread(fid,inf,'uint16')];
                    fclose(fid);
                end
            lastfullcycle = floor(length(DAT)/nch)*nch;
            DAT(lastfullcycle+1:end) = [];
            DAT = reshape(DAT,nch,[])';
            N = size(DAT,1); % samples
            S = fix(N/fs)-(winsz-1); % seconds
    
    % spectra
            PS = zeros(nch,S,length(freq));
            fprintf('computing spectra...');
            for ii = 1:S % compute in 'winsz'-sec windows in 1-sec steps
                fprintf('%02d%%',fix(ii/S*100));
                ind = round((ii-1)*fs+1:(ii+winsz-1)*fs);
                    for jj = 1:nch
                        PS(jj,ii,:) = periodogram(DAT(ind,jj),hamming(length(ind)),freq,fs); % TODO: remove mean before spectral analysis?
                    end
                if ii~=S, fprintf('\b\b\b'); end
            end
            fprintf('\n');
    
    % smoothing
            M = fix(S/60);
            PM = zeros(nch,size(PS,3),M);
            disp('smoothing spectra...');
                for ii = 1:M
                    ind = (ii-1)*60+1:ii*60;
                    for jj = 1:nch
                        PM(jj,:,ii) = trimmean(squeeze(PS(jj,ind,:)),10); % robust mean for every 60-sec (removes upper and lower 'winsz' seconds of each minute of data)
                    end
                end
            x = (1:M)/60;

            % save
            save([ddir subj filesep 'processed' filesep tank],'PM','x','freq');
        end

    % pre-plot
    ind = find(x>0.5 & x<21.5);
    x = x(ind);
    PM = PM(:,:,ind);
    for ich = 1:size(PM,1)
        PM(ich,:,:) = smooth2(squeeze(PM(ich,:,:)),2,1);
        PM(ich,:,:) = zscore(squeeze(PM(ich,:,:))')'; % zscore
    end

    % plot
    for ich = 1:size(PM,1)
        figure('Name',[tank ': ' num2str(ich)],'NumberTitle','off','Units','normalized','Position',[1/8 1/4 3/4 1/2],'Color','w');
        imagesc(x,1:length(freq),squeeze(PM(ich,:,:)));
        ytick = [1:10 20:10:100 200];
        yticknm = {'1';'';'';'';'';'';'';'';'';'10';'';'';'';'';'';'';'';'';'100';''};
        yind = zeros(size(ytick));
        for ii = 1:length(ytick)
            [~,yind(ii)] = min(abs(freq-ytick(ii)));
        end
        set(gca,'Box','off','YDir','normal','TickDir','out','YLim',[1 length(freq)],'XLim',[x(1) x(end)],'YTick',yind,'YTickLabel',yticknm);
        xlabel('time (h)'); ylabel('frequency (Hz)'); title([tank ': ' num2str(ich)],'Interpreter','none');
        set(gca,'CLim',[-3 6]);
        hc = colorbar; set(get(hc,'YLabel'),'String','power');
    end
end