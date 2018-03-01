function [C3_N2] = defineNREM(dat)

% Script to define NREM based on visual analysis for later recall by
% wamsley_spindles to detection spindles

NREMb = [8.5 9 10 11.6 16]; % hours for each begining segment of NREM
NREMe = [9 10 11 12.2 16.4]; % hours for each ending segment of NREM
fs = 1000; % sampling rate
timecon = 3600; % time conversion (hours from dataset)

for j = 1:length(NREMb);
    C3_N2{j} = dat(NREMb(j)*timecon*fs:NREMe(j)*timecon*fs);



end
