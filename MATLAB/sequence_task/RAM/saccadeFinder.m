function [Sindex, Ssize, Sdir] = saccadeFinder(x,y,fs)
% function [Sindex, Ssize, Sdir] = saccadeFinder(x,y,fs)
%   Find saccades from eyetracker data (x and y coordinates).
%   Sindex = N x 2 saccade index [start stop]
%   Ssize = N x 1 saccade magnitude (units of x and y)   
%   Sdir = N x 1 saccade direction (deg)
%
%   Andrew Murphy and Drew Richardson 10/2015

% parameters
sLength = 0.15; %sacade length in seconds, overestimates tend to be better
                %because the program tends to do a good job eliminating 
                %extraneous saccades
secondThresh = .1;   %second derivative trigger threshold in standard
                     %deviations (to avoid triggering off nois)
firstThresh = .3;   %first derivative trigger threshold
dns = 100; % downsampling factor
filt = 1; % filtering option
debug = 0; % debug option (0, 1, or 2)

% preprocessing
d = designfilt('lowpassfir','PassbandFrequency',.1,'StopbandFrequency',.2);
x = filtfilt(d,x); % filter
y = filtfilt(d,y);
if dns>1
    Fs = fs/100;
    x = downsample(x,dns); % downsample
    y = downsample(y,dns);
else
    Fs = fs;
end
sLengthFrms = round(sLength*Fs); % saccade length in frames
flength = Fs*0.02; % defining filter length in frames (20ms)
flength = 2.*round((flength+1)/2)-1;
if filt==1
    x = sgolayfilt(x,2,round(flength));
    x = medfilt1(x,round(flength));
    vx = diff(x); % first derivative: velocity
    y = sgolayfilt(y,2,round(flength));
    y = medfilt1(y,round(flength));
    vy = diff(y);
    vel = sqrt(vx.^2+vy.^2);
    vel = filtfilt(d,vel);
    acc = diff(vel); % second derivative: acceleration
    acc = filtfilt(d,acc);
    acc = sgolayfilt(acc,2,round(flength));
    acc = filtfilt(d,acc);
else
    vx = diff(x);
    vy = diff(y);
    vel = sqrt(vx.^2+vy.^2);
    acc = diff(vel);
end
N = length(acc);
x(N+1:end) = [];
y(N+1:end) = [];
vel(N+1:end) = [];

% artifact detection
ibad = (1:20)'; % indices of bad data (end effects of filtering and artifacts)
iart = find(abs(vel)>7*std(vel)); % artifact = huge changes in 1st derivative 
if ~isempty(iart) 
    swin = -sLengthFrms:sLengthFrms;
    iart = iart*ones(1,length(swin)) + ones(length(iart),1)*swin; % add window around art to bad indices list
    iart = unique(iart(:));
    iart(iart<1 | iart>N) = [];
    ibad = [ibad; iart];
end
if debug == 1
    dt = 5*Fs;
    t = (1:N)/Fs;
    for ii = 1:floor(N/dt)
        ind = round((ii-1)*dt+1):round(ii*dt);
        plot(t(ind),x(ind),'k','LineWidth',1); hold on;
        plot(t(ind),y(ind),'Color',[0.8 0.8 0.8],'LineWidth',1);
        iii = ismember(ind,ibad);
        if ~isempty(iii)
            plot(t(ind(iii)),x(ind(iii)),'ro');
            plot(t(ind(iii)),y(ind(iii)),'ro');
        end
        set(gca,'XLim',[t(ind(1)) t(ind(end))]);
        pause; close;
    end
end
igood = setdiff((1:N)',ibad);

% saccade detection
[ipk,ilc]=findpeaks(abs(acc));  %find peaks in the 2rd derivative
[vpk,vlc]=findpeaks(abs(vel));  %find peaks in the 1st derivative
iii = ismember(ilc,ibad);
ipk(iii) = []; ilc(iii) = [];   %remove peaks during artifacts
iiv = ismember(vlc,ibad);
vpk(iiv) = []; vlc(iiv) = [];
stdi = std(acc(igood));
stdv = std(vel(igood));
iii = find(abs(ipk)<secondThresh*stdi); % remove low amplitude peaks
ipk(iii) = []; ilc(iii) = [];
iiv = find(abs(vpk)<firstThresh*stdv);
vpk(iiv) = []; vlc(iiv) = [];
del = []; % remove velocity peaks that are too close to each other
for ii = 1:length(vlc)-1
    if vlc(ii+1)-vlc(ii)<sLengthFrms
        if vpk(ii)>vpk(ii+1)
            del = [del; ii+1];
        else
            del = [del; ii];
        end
    end
end
vpk(del) = []; vlc(del) = [];
del = []; % remove acceleration peaks that are too far from velocity peaks
for ii = 1:length(ilc)-1
    c = min(abs(vlc-ilc(ii)));
    if c > sLengthFrms/2
        del = [del; ii];
    end
end
ipk(del) = []; ilc(del) = [];
impvelN = zeros(length(vlc),1); % number of acceleration peaks for each velocity peak 
impvelLoc = zeros(length(ilc),1); % index of velocity peak closest to each acceleration peak
for ii = 1:length(ilc)
    [~,ind] = min(abs(vlc-ilc(ii)));
    impvelN(ind) = impvelN(ind) + 1;
    impvelLoc(ii) = ind;
end
del = [];
for ii = 1:length(vlc) % remove all but start and end acceleration peak for each velocity peak
    if impvelN(ii) > 2
        ind = find(impvelLoc==ii);
        ds = vlc(ii)-ilc(ind);
        [~,iii] = sort(ds);
        ind = ind(iii);
        del = [del; ind(2:end-1)]; 
    end
end
ipk(del) = []; ilc(del) = []; impvelLoc(del) = [];
Sindex = zeros(length(vlc),2);  % output: start/stop index
Ssize = zeros(length(vlc),1);   % output: saccade size
Sdir = zeros(length(vlc),1);    % output: saccade direction
del = [];
for ii = 1:length(vlc)
    ind = find(impvelLoc==ii);
    if length(ind)~=2, del = [del; ii]; continue; end
    Sindex(ii,:) = ilc(ind);
    Ssize(ii,:) = sqrt(diff(x(ilc(ind)))^2 + diff(y(ilc(ind)))^2);
    Sdir(ii,:) = atan2(diff(y(ilc(ind))),diff(x(ilc(ind))))*180/pi;
end
Sindex(del,:) = []; Ssize(del) = []; Sdir(del) = [];
if debug == 2
    dt = 5*Fs;
    t = (1:N)/Fs;
    figure('Name','saccade finder','NumberTitle','off','Units','normalized','Position',[1/8 1/8 3/4 3/4],'Color','w');
    for ii = 1:floor(N/dt)
        ind = round((ii-1)*dt+1):round(ii*dt);
        iii = find(Sindex(:,1)>ind(1) & Sindex(:,2)<ind(end));
        for jj = 1:length(iii)
            patch([t(Sindex(iii(jj),1)) t(Sindex(iii(jj),2)) t(Sindex(iii(jj),2)) t(Sindex(iii(jj),1))],[min(x(Sindex(iii(jj),1):Sindex(iii(jj),2)))*ones(1,2) max(x(Sindex(iii(jj),1):Sindex(iii(jj),2)))*ones(1,2)],'k','FaceColor','r','EdgeColor','none'); hold on;
            patch([t(Sindex(iii(jj),1)) t(Sindex(iii(jj),2)) t(Sindex(iii(jj),2)) t(Sindex(iii(jj),1))],[min(y(Sindex(iii(jj),1):Sindex(iii(jj),2)))*ones(1,2) max(y(Sindex(iii(jj),1):Sindex(iii(jj),2)))*ones(1,2)],'k','FaceColor','r','EdgeColor','none');
            text((t(Sindex(iii(jj),1))+t(Sindex(iii(jj),2)))/2,max(max(x(Sindex(iii(jj),1):Sindex(iii(jj),2))),max(y(Sindex(iii(jj),1):Sindex(iii(jj),2)))),[num2str(round(Ssize(iii(jj)))) ', ' num2str(round(Sdir(iii(jj))))],'HorizontalAlignment','center','VerticalAlignment','bottom');
        end
        plot(t(ind),x(ind),'k','LineWidth',1); hold on;
        plot(t(ind),y(ind),'Color',[0.8 0.8 0.8],'LineWidth',1); hold on;
        axis tight; set(gca,'XLim',[t(ind(1)) t(ind(end))]); xlabel('sec');
        pause; cla;
    end
end
if dns>1 % upsample saccade indices to get back to original sampling rate
    Sindex = (Sindex(:,1:2)-1)*dns+1;
end