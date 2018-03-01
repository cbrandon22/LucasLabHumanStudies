function plotSpindles(dat)
timecon = 3600;
winL = 8.5*timecon*fs;
winH = 9*timecon*fs;
% t = 1:length(dat);

figure
plot(t(winL:winH)/3600,zscore(EEGWave(winL:winH)));
 
hold on
plot(t(winL:winH)/3600,detection(winL:winH)*3);
hold on
% plot(t(winL:winH)/3600,zscore(dat(winL:winH)));

 

spindlecount = sum(detection(winL:winH)/1000)


end