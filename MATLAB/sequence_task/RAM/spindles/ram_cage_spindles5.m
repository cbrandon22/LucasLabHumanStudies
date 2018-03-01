% signalmean = threshold_wamsley(C3_N2,fs);
% detection = wamsley(C3,fs,signalmean);

    function signalmean = ram_cage_spindles
    
        ddir = '/Users/tnl/Desktop/'; % directory
        subj = 'fidel'; % subject
        tank = 'fidel_cage_111715'; % tank 
        fs = 1000; % sampling rate (Hz)
        timedo = 3600; % time correction 3600 = hours, 60 = minutes
        windL = 8.5*timedo; % data range
        windH = 14*timedo; % data range
        nch = 6; % number of channels
        dns = 50; % downsample factor
        %% Define parameters for the wavelet analysis
        fb = 13.5; % bandwidth parameter
        fc = 0.5; % wavelet center frequency
        scale = 3.7; % Freq = 10Hz, scale = 0.002: 250Hz, scale = 0.004: 125Hz
%         scal2frq(0.005,['cmor' num2str(fb) '-' num2str(fc)],1/fs) %Freq =
%         100Hz
        Ltotal = 0;

        winsz = 5; % window size for spectral estimates (s) - determines frequency resolution

        % load data
        % data
  % load data
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
            DAT = reshape(DAT,nch,[])';

            dat = DAT;
            N = size(dat,1); % samples
            S = fix(N/fs)-(winsz-1); % seconds

        if dns>1
            [b,a] = butter(2,round(fs/(2.5*dns))/(fs/2)); % antialiasing filter
            dat = filtfilt(b,a,dat);
        end
        dat = downsample(dat,dns);
        fs = round(fs/dns);
        t = (1:length(dat))/fs;
%         
%         T = N/fs; % sec


        % THRESHOLD_WAMSLEY Calculates the amplitude criteria for spindle detection
        % using the wamsley method.
        % Input is a cell of continuous segment of
        % EEG data recorded. The sampling frequency is the
        % final input.
        
            col = 3; % For Fidel Left CA3 = columns 3 and 4
            signal = DAT(windL:windH,col); 
            timee = 1:length(windL:windH);
            L = length(signal);
            %% Perform wavelet transformation
            figure(1);
            EEGWave = cwt(signal,scale,['cmor' num2str(fb) '-' num2str(fc)]);
            EEGData = real(EEGWave.^2);
            
            %% Take Moving Average
            EEGData = EEGData.^2;
            window = ones((fs/10),1)/(fs/10); % create 100ms window to convolve with
            EEGData2 = filter(window,1,EEGData); % take the moving average using the above window
            MA(Ltotal+1:Ltotal+L) = EEGData2;
            Ltotal = Ltotal+L;
            
            signalmean = mean(MA);
            
                    %% Determine amplitude threshold
        threshold = signalmean.*4.5; % defines the threshold
                %% Find Peaks in the MS Signal
        current_data=EEGData2;
        
        over=current_data>threshold; % Mark all points over threshold as '1'
        detection = zeros(length(current_data),1);
        detection(over) = 1;
        
        figure(1);
        
   
        plot(timee/timedo,zscore(signal));
        hold on
        plot(timee/timedo,detection*max(zscore(signal)));
        hold on
        plot(timee/timedo,zscore(EEGWave));
        
        [begins,ends] = find_spindles(detection);
        [detection,begins,ends] = maximum_duration(detection,begins,ends,3,fs);
        [detection,begins_03,ends_03] = minimum_duration(detection,begins,ends,0.3,fs);
        
        locs_03=(zeros(1,length(current_data)))';  % Create a vector of zeros the length of the MS signal
        for i=1:((length(current_data))-(fs*0.01));  % for the length of the signal, if the sum of 30 concurrent points = Fs*0.3, mark a spindle
            if sum(over(i:(i+((fs*0.3)-1))))==(fs*0.3);
                locs_03(i,1)=1;
            end
        end
        
        spin_03=zeros((length(locs_03)),1);  % only mark a spindle in vector 'spin' at the end of a 300ms duration peak
        for i=1:length(locs_03);
            if locs_03(i,1)==1 && locs_03(i+1,1)==0;
                spin_03(i,1)=1;
            end
        end
        
        figure
        
        figure (3);
        plot(timee/timedo,EEGWave);
        hold on
        plot(timee/timedo,spin_03*max(EEGWave));
        
%         for i=201:length(spin_03);  % for every spindle marked in 'spin', delete the spindle if there is also a spindle within the second preceeding it
%             if spin_03(i,1)==1 && sum(spin_03((i-fs):(i-1)))>0;
%                 spin_03(i,1)=0;
%                 idx = find(i>=begins_03 & i<=ends_03);
%                 if isempty(idx) == 0
%                     detection(begins_03(idx):ends_03(idx)) = 0;
%                 else
%                     error('Did not find spindle beginning and ending around a spin point')
%                 end
%             end
%         end
        
%           save([ddir subj filesep 'processed' filesep tank '_vt2'],'dat','EEGWave','EEGFreq','EEGsgram','EEGData2','t','N','S','fs','detection','begins','ends');
    end

  
    
