function [begins_03,ends_03,detection,signal] = rippledetect(dat);

% Bandpass filter LFP signal 150-200Hz
%         col = 1; % For Fidel Left CA3 = columns 3 and 4, ditka is always 1
        fs = 10000;
        ripple_duration = [0.025 0.15]; %ripple duration/width = 0.02 - 0.15 seconds or 20-150ms.
        thresh = 4; %amplitude threshold floor
        noise1 = 7; % amplitude threshold ceiling to minimize large amplitude noise
        datastart = 1;
        dataend = 8.4;  
        timedo = 3600*fs; % time correction 3600 = hours, 60 = minutes
        windL = datastart; % data range
        windH = dataend*timedo; % data range
        Ltotal = 0;
        if datastart == 1
            timecorrection = 1; % convert to allow hours entry above, use 1 if datastart = 1;
        else 
            timecorrection = 3600*fs; % convert to allow hours entry above, use 1 if datastart = 1;
        end

        data = dat(datastart*timecorrection:dataend*3600*fs);
        [b,a] = butter(2,[150 200]/(fs/2)); % Ripple band filter 150-200 Hz
        data = filtfilt(b,a,data);

% Calculate RMS, moving 10ms window, Ripple detection 4 SD above
% RMS signal (see threshold and noise levels below).

            %% Take Moving Average
            signal = data; 
            L = length(signal);
            EEGData = signal.^2;
            window = ones((fs/100),1)/(fs/100); % create 10ms window? to average (assumes 1000Hz sampling frequency
            EEGData2 = filter(window,1,EEGData); % take the moving average using the above window
            MA(Ltotal+1:Ltotal+L) = EEGData2;
            Ltotal = Ltotal+L;
            
            signalmean = mean(MA);
                    %% Determine amplitude threshold
        threshold = signalmean + (thresh*std(signal)); % defines the threshold
        noise = signalmean + (noise1*std(signal)); % defines noise level
                %% Find Peaks in the MS Signal
        current_data=EEGData2;
        
        over=current_data>threshold & current_data<noise; % Mark all points over threshold as '1'
        detection = zeros(length(current_data),1);
        detection(over) = 1;
        %% Time/Window ripple determiniation. ripple width = 0.02 - 0.15 seconds or 20-150ms.
        [begins,ends] = find_spindles(detection);
        [detection,begins,ends] = maximum_duration(detection,begins,ends,ripple_duration(2),fs);
        [detection,begins_03,ends_03] = minimum_duration(detection,begins,ends,ripple_duration(1),fs);
        
        figure('Name','Ripples');
        
   
        plot((windL:1:windH)/timedo,signal);
        hold on
        plot((windL:1:windH)/timedo,detection*max(signal));
        xlabel({'Time (h)'});
        ylabel({'Amplitude'});
        title({'Ripples Detected'});
        
%                 locs_03=(zeros(1,length(current_data)))';  % Create a vector of zeros the length of the MS signal
%         for i=1:((length(current_data))-(fs*0.01));  % for the length of the signal, if the sum of 10 concurrent points = Fs*0.01, mark a spindle
%             if sum(over(i:(i+((fs*0.01)-1))))==(fs*0.01);
%                 locs_03(i,1)=1;
%             end
%         end
%         
%         spin_03=zeros((length(locs_03)),1);  % only mark a ripple in vector 'spin' at the end of a 300ms duration peak
%         for i=1:length(locs_03);
%             if locs_03(i,1)==1 && locs_03(i+1,1)==0;
%                 spin_03(i,1)=1;
%             end
%         end
%         plot((windL:1:windH)/timedo,spin_03*max(signal));
        
end