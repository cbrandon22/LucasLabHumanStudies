subj = 'HUP121_i';
timeWin = 500;
timeStep = 100;
freQ = logspace(log10(2), log10(200),50);
freqBins = logspace(log10(2), log10(200),50);

cd ~/../../Volumes/Lucas_ECoG_Data/Consciousness/data/;
cd (subj);
cd lfp.noreref;
filename = 'HUP121_i_01Jul16_0822.001';
[dat,samplerate] = loadeData(filename);
periodogram(dat);
figure;
cwt(dat,freQ,'morl');
% [~,pow] = multiphasevec2(freQ,dat',samplerate,7);
% % log before mean
% pow = log10(pow);
% % create time bins
% timeBins = timeWin:timeStep:size(pow,2);
% 
% %Initialize power mat 
% binPowMat = nan([length(freqBins) length(timeBins) size(pow,3)]);
% % loop through freqBins
% for f = 1:length(freqBins)
%     %identify freq of interest
%     [~,fIndEnd] = min(abs(freQ-freqBins(f)));
%     if f == 1
%         fIndBegin = 1;
%     else
%         [~,fIndBegin] = min(abs(freQ-freqBins(f-1)));
%         fIndBegin = fIndBegin+1; %b/c the closest value to bin edge was already
%                                  %used in the last iteration
%     end
%     fInd = fIndBegin:fIndEnd;
%     %loop through timeBins
%     for t = 1:length(timeBins)
%         
%         %identify time windows
%         tInd = (timeBins(t)-timeWin+1):timeBins(t);
%         
%         %calculate mean power
%         binPowMat(f,t,:) = nanmean(nanmean(pow(fInd,tInd,:),1),2);
%     end
%     
% end
% timeMinutes = (length(dat)/samplerate)/60;
% imagesc(0:timeMinutes,freqBins,binPowMat);
% keyboard
