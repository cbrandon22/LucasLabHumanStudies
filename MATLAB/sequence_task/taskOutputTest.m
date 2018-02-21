timestamp = xlsread('tsTestLog.xlsx','A1:A1000');
[~,message] = xlsread('tsTestLog.xlsx','B1:B1000');

%responseTones = []; %RESP,tone
%recallTones = [];
%trialStart = [];
cashSound = [];
for i=1:length(message)
    if strcmp(message{i},'cashtime')
        cashSound = [cashSound;timestamp(i),timestamp(i+1)];
    end
%     if strcmp(message{i},'RESP')
%         responseTones = [responseTones;timestamp(i),timestamp(i-1)];
%     elseif strcmp(message{i},'RECALL')
%         recallTones = [recallTones;timestamp(i),timestamp(i+1)];
%     elseif strcmp(message{i},'TRIAL_START')
%         trialStart = [trialStart; timestamp(i),timestamp(i+1),timestamp(i-1),timestamp(i-4)];
%     elseif strcmp(message{i},'RECALL_END') && strcmp(message{i-2},'soundtime')
%         cashSound = [cashSound;timestamp(i),timestamp(i-1)];
%     end
end

% resp2sound = responseTones(:,1)-responseTones(:,2); % time from RESP to actual tone
% recall2Tone = recallTones(:,2)-recallTones(:,1); % time from recall to actual tone
% trialstartdiff = trialStart(:,2)-trialStart(:,1);
% trialstarttest = trialStart(:,1)-trialStart(:,4);
cashSoundLength = cashSound(:,2)-cashSound(:,1);
