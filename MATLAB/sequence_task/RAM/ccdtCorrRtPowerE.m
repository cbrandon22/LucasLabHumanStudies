function ccdtCorrRtPowerE
% function ccdtCorrRtPowerE
%   Neural correlates of CCDT: neural power in a chosen frequency band and
%   time window as a function of movement time on error trials.
%
%   DR 03/2016

% parameters
ddir = '/Users/tnl/Desktop/'; % directory
subj = 'fidel'; % subject
dday = [110415; 110615; 111115; 111315; 111615]; % day(s)
ch = 1:16; % channels
band = [5 10]; % frequency band (Hz)
twin = [1400 1700]; % cue-aligned window (ms)
dtbins = 850:100:1450;  

% compile data
cd([ddir subj '/processed/ccdtCorr/']);
M = length(dtbins)-1;
Pbin = zeros(length(ch),M);
for ich = 1:length(ch)
    Pa = []; DTa = [];
    for iday = 1:length(dday)
        fls = dir([subj '_*_' num2str(dday(iday)) '_*_' num2str(ch(ich)) '.mat']);
        if length(fls)~=1, disp('file not found'); return; end
        load(fls(1).name);
        DTa = [DTa; DTe];
        Pa = [Pa; abs(Pcce)]; % amplitude
    end
    
    % window on cue event and compute average power
    if ich==1
        fs = 1000/mean(diff(t)); % samp/s
        swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
        ifreq = find(freq>=band(1)&freq<=band(2));
    end
    Pavg = NaN*ones(size(Pa,1),1);
    for ii = 1:size(Pa,1)
        [~,ind] = min(abs(t+DTa(ii)));
        try
            c = squeeze(Pa(ii,ifreq,swin+ind));
            Pavg(ii) = trimmean(c(:),5);
        end
    end
    
    % binned power
    for ii = 1:M
        ind = find(DTa>=dtbins(ii) & DTa<dtbins(ii+1));
        Pbin(ich,ii) = trimmean(Pavg(ind),5);
    end
end

% plot
figure('Name','ccdtCorrRTPower','NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
plot(dtbins(1:end-1)+mean(diff(dtbins))/2,Pbin-Pbin(:,1)*ones(1,M));
set(gca,'Box','off');
xlabel('movement time (ms)'); ylabel([num2str(band(1)) '-' num2str(band(2)) ' Hz power']);
