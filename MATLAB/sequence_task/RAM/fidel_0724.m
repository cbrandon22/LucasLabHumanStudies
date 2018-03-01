function fidel_0724
% function fidel_0724

% parameters
ddir = 'G:\'; % directory
subj = 'fidel'; % subject
tank = 'fidel_VPLT_DT1_072415'; % tank
block = 4; % block
chrec = 8; % rec channel(s)

% main
cd([ddir subj '\' tank '\Block-' num2str(block)]);
fid = fopen([tank '_Block-' num2str(block) '_xRAW_ch25.sev'],'r');
fseek(fid,40,'bof');
dat = fread(fid,inf,'single');
fclose(fid);
fs = 24414.06;
itrain = find(dat(1:end-1)<0.5 & dat(2:end)>0.5); % trigger for each pulse train
itrain = itrain + round(0.5/1000*fs);
iinfo = [itrain(end-19), 2, -1000, 8500; % index, freq, twin(1), twin(2)
        itrain(52), 4, -500 4250;
        itrain(end-17), 6, -333, 2833;
        itrain(end-15), 8, -250, 2125;
        itrain(end-12), 10, -500, 3000;
        itrain(end-9), 12, -500, 3000;
        itrain(end-6), 14, -500, 3000;
        itrain(end-5), 16, -500, 3000;
        itrain(end-3), 18, -500, 3000;
        itrain(end), 20, -500, 3000];
% iinfo = [itrain(end-19), 2, -50+7500, 150+7500; % index, freq, twin(1), twin(2)
%         itrain(52), 4, -50+3750, 150+3750;
%         itrain(end-17), 6, -50+2500, 150+2500;
%         itrain(end-15), 8, -50+1875, 150+1875;
%         itrain(end-12), 10, -50+1900, 150+1900;
%         itrain(end-9), 12, -50+2418, 150+2418;
%         itrain(end-6), 14, -50+2072, 150+2072;
%         itrain(end-5), 16, -50+1813, 150+1813;
%         itrain(end-3), 18, -50+1611, 150+1611;
%         itrain(end), 20, -50+1451, 150+1451];
fid = fopen([tank '_Block-' num2str(block) '_xRAW_ch' num2str(chrec) '.sev'],'r');
fseek(fid,40,'bof');
dat = fread(fid,inf,'single'); % uV
fclose(fid);
for ii = 1:size(iinfo,1)
    figure('Name',[num2str(iinfo(ii,2)) ' Hz'],'NumberTitle','off','Units','normalized','Position',[1/4 1/4 1/2 1/2],'Color','w');
    swin = round(iinfo(ii,3)/1000*fs):round(iinfo(ii,4)/1000*fs);
    indwin = iinfo(ii,1) + swin;
    t = swin/fs*1000;
    datwin = dat(indwin);
    datwin = datwin - median(datwin); % remove dc offset
    plot(t,datwin,'k');
    set(gca,'Box','off','YLim',[-1000 1000],'XLim',[t(1) t(end)]);
    xlabel('ms'); ylabel('\muV'); title([num2str(iinfo(ii,2)) ' Hz']);
end