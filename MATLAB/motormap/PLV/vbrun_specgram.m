function [Sg, t, fg] = vbrun_specgram(data, fs)

% Variables
dsf = 1; %should match the corresponding extract function dsf
fsd = fs/dsf; %should match the corresponding extract function fsd

params.Fs = fsd;
params.fpass = [1 120];
params.err = [2 .05];
params.tapers = [7 13];
params.trialave=1; % if averaging over multiple channels
movingwin = [.5 .25];


[Sg,t,fg]=mtspecgramc(data,movingwin,params);

% also compute fft alone
Y=fft(data);
L=length(data);
P2=abs(Y/L);
P1=P2(1:L/2+1);
P1(2:end-1)=2*P1(2:end-1);
f=fs*(0:(L/2))/L;

% plot
figure
plot(f,P1);
xlim([0 120])

figure
pcolor(t,fg,log(Sg'))
shading flat
colorbar