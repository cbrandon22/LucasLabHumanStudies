dxls = 'F:\ditka\ditka.xlsx'; % log
tsurg = [2015 03 03 00 00 00]; % surgery date

sess = xlsread(dxls,'IMP','A5:A14');
HCm = xlsread(dxls,'IMP','J5:J14');
HCv = xlsread(dxls,'IMP','K5:K14');
ECm = xlsread(dxls,'IMP','J18:J27');
ECv = xlsread(dxls,'IMP','K18:K27');
MSm = xlsread(dxls,'IMP','J31:J40');
MSv = xlsread(dxls,'IMP','K31:K40');

sess = num2str(sess);
sess = [str2double(cellstr(sess(:,1:4))), str2double(cellstr(sess(:,5:6))), str2double(cellstr(sess(:,7:8))), zeros(size(sess,1),3)];
x = etime(sess,ones(size(sess,1),1)*tsurg)/3600/24; % days since surgery
figure('Name','ditka impedance','NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');    
errorbar(x,HCm,HCv,'b'); hold on;
errorbar(x,ECm,ECv,'r');
errorbar(x,MSm,MSv,'g');
axis square; set(gca,'Box','off');
xlabel('days post implant'); ylabel('impedance at 1000 Hz (k\Omega)')
legend({'hippocampus','entorhinal cortex','medial septum'});
legend('boxoff');
