function [stplv] = stPLV(eegData, srate, filtSpec, twin)
% Computes the Single-Trial Phase Locking Value (PLV) for an EEG dataset.
%
% Input parameters:
%   eegData is a 3D matrix numChannels x numTimePoints x numTrials
%   srate is the sampling rate of the EEG data
%   filtSpec is the filter specification to filter the EEG signal in the
%     desired frequency band of interest. It is a structure with two
%     fields, order and range.
%       Range specifies the limits of the frequency
%     band, for example, put filtSpec.range = [35 45] for gamma band.
%       Specify the order of the FIR filter in filtSpec.order. A useful
%     rule of thumb can be to include about 4 to 5 cycles of the desired
%     signal. For example, filtSpec.order = 50 for eeg data sampled at
%     500 Hz corresponds to 100 ms and contains ~4 cycles of gamma band
%     (40 Hz).
%   twin - sliding time window to average instantaneous phase difference (s)
%
% Output parameters:
%   stplv is a structure with:
%   plv is a 4D matrix -
%     numTimePoints x numChannels x numChannels x numTrials
%   et is the time estimate per trial
%   trial is the trial number
%
%--------------------------------------------------------------------------
%
% NOTE:
% As you have probably noticed in the plot from the above example, the PLV
% between two random signals is spuriously large in the first 100 ms. While
% using FIR filtering and/or hilbert transform, it is good practice to
% discard both ends of the signal (same number of samples as the order of
% the FIR filter, or more).
%
% Also note that in order to extract the PLV between channels 17 and 20,
% use plv(:, 17, 20, :) and NOT plv(:, 20, 17, :). The smaller channel
% number is to be used first.
%--------------------------------------------------------------------------
%
% Reference:
%   Lachaux, J P, E Rodriguez, J Martinerie, and F J Varela.
%   �Measuring phase synchrony in brain signals.�
%   Human brain mapping 8, no. 4 (January 1999): 194-208.
%   http://www.ncbi.nlm.nih.gov/pubmed/10619414.
%
%--------------------------------------------------------------------------
% Written by:
% Vivek P. Buch, MD
% Neurosurgery Resident
% University of Pennsylvania

% Modified from PLV.m written by:
% Praneeth Namburi
% Cognitive Neuroscience Lab, DUKE-NUS
% 01 Dec 2009
%

numChannels = size(eegData, 1);
numTrials = size(eegData, 3);
stplv = struct;

% disp('Filtering data...');
filtPts = fir1(filtSpec.order, 2/srate*filtSpec.range);
filteredData = filter(filtPts, 1, eegData, [], 2);

% disp(['Calculating PLV for ' mat2str(sum(dataSelectArr, 1)) ' trials...']);
for channelCount = 1:numChannels
    filteredData(channelCount, :, :) = angle(hilbert(squeeze(filteredData(channelCount, :, :))));
end
swin = floor(srate*twin);
numTimeWindows = floor(size(eegData,2)/swin);
tWindows = zeros(numTimeWindows,swin);
for x = 1:numTimeWindows
    tWindows(x,:) = ((x-1)*swin)+1:((x-1)*swin)+swin;
end
for trialCount = 1:numTrials
    tic;
    fprintf(['trial ' num2str(trialCount) '/' num2str(numTrials) '\n'])
    tmpplv = zeros(size(filteredData,2),numChannels, numChannels);
    for channelCount = 1:numChannels-1
        channelData = squeeze(filteredData(channelCount, :, trialCount));
        for compareChannelCount = channelCount+1:numChannels
            compareChannelData = squeeze(filteredData(compareChannelCount, :, trialCount));
            for tCount = 1:numTimeWindows
                currTwin = tWindows(tCount,:);
                tmpplv(currTwin, channelCount, compareChannelCount) = abs(sum(exp(1i*(channelData(currTwin) - compareChannelData(currTwin))),2))/length(currTwin);
            end
        end
    end   
  stplv(trialCount).plv = tmpplv;
  stplv(trialCount).trial = trialCount;
  stplv(trialCount).et = toc;
  clear tmpplv
end
return;