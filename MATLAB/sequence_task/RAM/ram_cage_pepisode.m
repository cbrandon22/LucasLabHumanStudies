function [ep,t] = ram__cage_pepisode(varargin)
% function [ep,t] = ram_pepisode(varargin)
%   Oscillatory episode detection using Pepisode method (Caplan et al
%   2001).
%
%   DR 03/2015

% parameters
    ddir = '/Users/tnl/Desktop/'; % directory
    subj = 'ditka'; % subject
    tank = 'ditka_cage_040615'; % tank 
    fs = 1000; % sampling rate (Hz)
%     block = 1; % block
    nch = 1; % number of channels
    dns = 50; % downsample factor
    wavl = 'morl'; % wavelet
    fev = 11:1:16; % evaluation frequencies
    tcyc = 3; % duration threshold (cycles)
if nargin
    v2struct(varargin{1});
end

% data
 cd([ddir subj filesep tank]);
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
    dat = reshape(DAT,nch,[])';
    
if dns>1
    [b,a] = butter(2,round(fs/(2.5*dns))/(fs/2)); % antialiasing filter
    dat = filtfilt(b,a,dat);
end
dat = downsample(dat,dns);
fs = round(fs/dns);
t = (1:length(dat))/fs;

% wavelet transform
freq = logspace(log10(3),log10(100),50)'; % spectrum frequencies (Hz) %%%SENSITIVE
scales = centfrq(wavl)./(freq*(1/fs));
P = cwt(dat,scales,wavl);
P = abs(P); % power spectra

% power criterion
Pm = median(P,2); % mean power spectrum %%%MEDIAN
b = regress(log10(Pm),[ones(size(freq)) log10(freq)]); % colored noise fit [P(f) = log10(b(1))*f^(-b(1))]
ep = zeros(length(fev),length(Pm));
for ii = 1:length(fev)
    [~,ifev] = min(abs(freq-fev(ii)));
    mfev = 10^(b(1)+b(2)*log10(freq(ifev))); % mean power at evaluation frequency
    Pt = ncx2inv(.95,2,mfev); % power threshold (uV)

    % figure: power threshold
    if ii == 1
        figure
        h(1) = plot(freq,Pm,'k','LineWidth',2);
        set(gca,'Box','off','XScale','log','YScale','log'); hold on;
        h(2) = line(get(gca,'XLim'),10.^(b(1)+b(2)*log10(get(gca,'XLim'))),'Color',[0.8 0.8 0.8],'LineWidth',2);
        ylabel('power'); xlabel('Hz'); axis tight;
        legend(h,'mean spectrum','background spectrum'); legend('boxoff');
    end
    line(freq(ifev)*ones(1,2),[mfev Pt],'Color','r');
    plot(freq(ifev),Pt,'ro','MarkerFaceColor','r');
    
    % duration criterion
    Dt = fix(tcyc/freq(ifev)*fs); % duration threshold (samples)
    Dt = Dt + mod(Dt-1,2); % ensure it is odd number of samples so that length(Dtx) = Dt
    Dtx = -floor(Dt/2):floor(Dt/2); % samples before and after center of an episode
    
    % episodes
    cep = filter(ones(Dt,1),1,double(P(ifev,:)>Pt))==Dt; % moving sum filter with Dt/2 sample delay
    iep = find(cep)'; iep = iep-fix(Dt/2); % correct centers for delay
    iepmat = (iep*ones(1,Dt) + ones(length(iep),1)*Dtx)'; 
    ep(ii,iepmat(:)) = 1; % make all Dt samples = 1, not just center
end
