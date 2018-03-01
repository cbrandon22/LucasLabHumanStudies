function [pks,iRip,detection,signal] = rippledetectFP(dat);

% Bandpass filter LFP signal 150-200Hz
%         col = 1; % For Fidel Left CA3 = columns 3 and 4, ditka is always 1
        fs = 1000;
        ripple_duration = [0.025 0.15]; %ripple duration/width = 0.025 - 0.15 seconds or 25-150ms.
        th = 4; %amplitude threshold floor
        datastart = 8;
        dataend = 17;  
        timedo = 3600*fs; % time correction 3600 = hours, 60 = minutes
        windL = datastart*timedo; % data range
        windH = dataend*timedo; % data range

        if datastart == 1
            timecorrection = 1; % convert to allow hours entry above, use 1 if datastart = 1;
        else 
            timecorrection = 3600*fs; % convert to allow hours entry above, use 1 if datastart = 1;
        end

        data = dat(datastart*timecorrection:dataend*3600*fs);
        
        % ripple detection
        widthMin = ripple_duration(1);
%         widthMax = ripple_duration(2);
        swidMin = round(widthMin/1000*fs);
%         swidMax = round(widthMax/1000*fs);
        [bB,aB] = butter(2,[150 250]/(fs/2)); %ripple bandpass filter
        [bA,aA] = butter(2,20/(fs/2)); % lowpass averaging filter
        signal = filtfilt(bB,aB,data); % ripple bandpass filtered signal
        Rrms = filtfilt(bA,aA,abs(signal)); % RMS of ripple-band activity
        
        thresh = mean(Rrms)+th*std(Rrms);
        [pks,iRip] = findpeaks(Rrms,'MinPeakHeight',thresh,'MinPeakWidth',swidMin); % ,'MaxPeakWidth',swidMax
        detection = zeros(length(signal),1);
        detection(iRip) = 1;

        %% Time/Window ripple determiniation. ripple width = 0.02 - 0.15 seconds or 20-150ms.
%         [begins,ends] = find_spindles(iRip);
     
        figure('Name','RipplesFP');
           
        plot((windL:1:windH)/timedo,signal);
        hold on
        plot((windL:1:windH)/timedo,detection*max(signal));
        xlabel({'Time (h)'});
        ylabel({'Amplitude'});
        title({'Ripples Detected'});

        
end