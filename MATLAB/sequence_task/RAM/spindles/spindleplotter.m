function spindleplotter
% [spindles_start_end, detected_spindles, coefs] = spindle_estimation_FHN2015(dat(3600*1000:2*3600*1000), 1000);

fs = 1000;
Spind = 150;
win1 = 500;
win2 = 500;
Volts1 = dat(spindles_start_end(Spind,1)-win1:spindles_start_end(Spind,2)+win2);
% VoltsBPF = datBPFSpindles(spindles_start_end(Spind,1):spindles_start_end(Spind,2));
VoltsZ = zscore(Volts1);
Time1 = t(spindles_start_end(Spind,1)-win1:spindles_start_end(Spind,2)+win2);
Time2 = 0-win1:spindles_start_end(Spind,2)+win2-spindles_start_end(Spind,1);
% Time3 = Time1-spindles_start_end(Spind,2)+win2;
figure('Name',['EEG Spindle #' Spind '(D040615)'])

plot(Time2/1e3,VoltsZ);
hold on
plot(Time1/1e6,zscore(VoltsBPF));

        bpf = [11 16]; % bandpass filter freq (Hz) - 11 16spindle frequencies
%         bpf2 = [140 200]; % bandpass filter freq (Hz) - ripple frequencies
        fs = 1000; % sampling rate (Hz)
    [b,a] = butter(2,bpf/(fs/2));
%     [b2,a2] = butter(2,bpf2/(fs/2));
    datBPFSpindles = filtfilt(b,a,Volts1);
%     datBPFRipples = filtfilt(b2,a2,DAT);
figure('Name','Bandpass Filtered Spindle')

plot(Time2/1e3,datBPFSpindles);

% wname = 'morl'; % Morlet wavelet 'morl'
% scales = 2:0.1:15; % scales corresponding to pseudo-frequencies
% VoltsCWT = cwt(Volts1,scales,wname); % Apply CWT 
% % VoltsCWT = coefs;
% VoltsCWTA = mean(VoltsCWT);
        fb = 13.5;
        fc = 0.5;
        scale = 3.7;
EEGWave = cwt(Volts1,scale,['cmor' num2str(fb) '-' num2str(fc)]);




hold on
plot(Time1/3600/fs,EEGWave);


end