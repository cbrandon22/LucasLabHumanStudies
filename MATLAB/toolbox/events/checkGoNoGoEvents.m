% Script specifically written to check the events files of 2 subjects for
% Xenon study using GoNogo task. Concatenates events across sessions

subjDir = '/Users/tnl/Desktop/GoNoGo Subject 3 Output';
i=0;
while isdir(fullfile(subjDir,sprintf('session_%s',num2str(i))))
    load(fullfile(subjDir,sprintf('session_%s',num2str(i)),'events.mat'),'events')
    if exist('goodEv','var')
        goodEv = [goodEv events];
    else
        goodEv = events;
    end
    i=i+1;
end

subjDir = '/Users/tnl/Desktop/xenon4';
i=0;
while isdir(fullfile(subjDir,sprintf('session_%s',num2str(i))))
    load(fullfile(subjDir,sprintf('session_%s',num2str(i)),'events.mat'),'events')
    if exist('badEv','var')
        badEv = [badEv events];
    else
        badEv = events;
    end
    i=i+1;
end

%Get start of each session
gSessStart = find(strcmp({goodEv.type},'B'));
bSessStart = find(strcmp({badEv.type},'B'));

%check session lengths
gSessLen = diff(gSessStart);
bSessLen = diff(bSessStart);

%Trim short session from badEv, confirm lengths again
origBadEv = badEv;
badEv = badEv(bSessStart(2):end);
bSessStart = find(strcmp({badEv.type},'B'));
bSessLen = diff(bSessStart);

%Get Response Starts
gRespStart = find(strcmp({goodEv.type},'RESPONSE_START'));
bRespStart = find(strcmp({badEv.type},'RESPONSE_START'));
gstart2first = cell2mat({goodEv(gRespStart+1).mstime})-cell2mat({goodEv(gRespStart).mstime});
bstart2first = cell2mat({badEv(bRespStart+1).mstime})-cell2mat({badEv(bRespStart).mstime});
plot([1:1:length(gstart2first)],gstart2first,'o',[1:1:length(bstart2first)],bstart2first,'o')
