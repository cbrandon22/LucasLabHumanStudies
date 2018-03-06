%% Filter Comparisons
%
% Currently getting different results between idealfilt and filtfilt
% functions. Want to use both on some simulated data to see what the
% difference is
%
% @author JStiso

%% Test Low and High Pass

% set filter interval
lpf = 40;
hpf = 1;
srate = 1000;

% initialize ecg data
wform = ecg(500);
x = wform' + 0.25*randn(500,1);

% make filter
[a, b] = butter(4, [hpf/(srate/2), lpf/(srate/2)]);
% filtfilt
y = filtfilt(a, b, x);
% phase distorting filter
y1 = filter(a, b, x);
% ideal filt
m = mean(x);
x_m = x - m;
x_ts = timeseries(x_m, (1:numel(x))./srate);
y_if = idealfilter(x_ts, [hpf lpf], 'pass');
y_if.Data = y_if.Data + m;

subplot(2,1,1)
plot([y y_if.Data y1])
title('Filtered Waveforms')
legend('Butterworth', 'Idealfilter', 'Conventional Filtering')

subplot(2,1,2)
plot(x)
title('Original Data')

figure;
clf;
subplot(3, 1, 1)
spectopo(y', 0, 1000)
title('Butterworth')
subplot(3, 1, 2)
spectopo(y_if.data', 0, 1000)
title('Idealfilter')
subplot(3, 1, 3)
spectopo(y', 0, 1000)
title('Conventional Filter')

%% Test Notch Filter

% set filter interval
low_cut = 30;
high_cut = 100;
srate = 1000;

% initialize ecg data
wform = ecg(500);
x = wform' + 0.25*randn(500,1);

% make filter
[a, b] = butter(4, [low_cut/(srate/2), high_cut/(srate/2)], 'stop');
% filtfilt
y = filtfilt(a, b, x);
% phase distorting filter
y1 = filter(a, b, x);
% ideal filt
m = mean(x);
x_m = x - m;
x_ts = timeseries(x_m, (1:numel(x))./srate);
y_if = idealfilter(x_ts, [low_cut, high_cut], 'notch');
y_if.Data = y_if.Data + m;

subplot(2,1,1)
plot([y y_if.Data y1])
title('Filtered Waveforms')
legend('Filtfilt', 'Idealfilter', 'Conventional Filtering')

subplot(2,1,2)
plot(x)
title('Original Data')

figure;
clf;
subplot(3, 1, 1)
spectopo(y', 0, 1000)
title('Butterworth')
subplot(3, 1, 2)
spectopo(y_if.data', 0, 1000)
title('Idealfilter')
subplot(3, 1, 3)
spectopo(y', 0, 1000)
title('Conventional Filter')
