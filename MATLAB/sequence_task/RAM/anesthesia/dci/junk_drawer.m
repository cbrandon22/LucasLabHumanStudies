
%% Simulate Junk Drawer Effect
% This script is meant to explore the "junk drawer effect" observed in my
% human ECoG data. Essentially it seems like the model needs access to high
% frequencies in order to categorize unstable (decorrelated across space)
% modes. Without it, it is only capable of finding stable modes.
%
% The plan here is to simulate an AR process with eigenvalues that follow
% some power law. Then take the simulated data, filter it at different low
% pass cut offs, and see for which frequencies it is possible to recover
% the initial eigenmode spectrum
%
% @author JStiso April 2017

%% Define constants

n_elecs = 128;
srate = 1000;
% initialize frequencies
freq = zeros(1, n_elecs);
% add negatives to get eigenvalues in conjugate pairs
freq_pos = logspace(log10(3), log10(500), n_elecs/2); % in hz
freq(1:2:n_elecs - 1) = freq_pos;
freq(2:2:n_elecs) = -1.*freq_pos;
pow_exp = -1.1;
dampr = -1./(abs(freq).^pow_exp); % in 1/s

% should look like a power law
clf
scatter(dampr, freq);

%% Make eigenvalues
% lambda = p*exp(i*phi)
% dampr = log2(p)*srate
% freq = phi*srate/2pi

p = (dampr./srate).^2;
phi = (freq.*2.*pi)/srate;

% lambda = a + bi
% lambda = pcos(phi) + isin(phi)

lambda = (p.*(cos(phi) + 1i.*sin(phi)))';

% values should be less than one
clf
plot(abs(lambda))


% make A
A = diag(lambda);

%% Simulate AR1 process

% constants
time = 2000;
x = zeros(n_elecs,time);
x(:,1) = rand(1,n_elecs);

for t = 1:time-1
    x(:,t+1) = A*x(:,t) + rand;
end
x = abs(x);

plot(x(1,1:500))
%eegplot(abs(x), 'srate', srate)

%% Filter data, and try to recover eigenvalues

% window size
win = 500;
filters = 100:50:500;

x_filt = zeros(size(x));
for f = 1:numel(filters)
    
    %use ideal filter
    % for every electrode make x a timeseries object
    for e = 1:size(x,1)
        curr_x = x(e,:);
        % demean for idealfilter
        m = mean(curr_x,2);
        curr_x = curr_x - m;
        x_ts = timeseries(curr_x,(1:size(curr_x,2))./srate);
        % high and low pass filters
        interval = [0, filters(f)];
        x_ts = idealfilter(x_ts, interval, 'pass');
        % make array again
        tmp = x_ts.data;
        % add mean back
        tmp = tmp + m;
        x_filt(e,:) = tmp;
    end
    
    %initialize
    eigvect = zeros(size(x_filt,1), size(x_filt,1), size(x_filt,2) - win);
    eigval = zeros(size(x_filt,1), size(x_filt,2) - win);
    parfor i = 1:size(x_filt,2) - win
        fprintf('\n...%d', i)
        x_curr = x_filt(:,i:i + win);
        [~,A,~,~,~,~] = arfit(x_curr',1,1,'sbc'); % data_n should be a time chunk;
        [vect, val] = eig(A);
        eigvect(:,:,i) = vect;
        eigval(:,i) = diag(val);
    end
    
    
    %vectorize, get abs, and get dist
    ref_em1 = reshape(eigval(:,1:time-win), [], 1);
    [n1, x1] = hist(abs(ref_em1),100);
    
    %get dist of true evs
    [n2, x2] = hist(abs(lambda), 100);
    
    % normalize
    n1 = n1/max(n1);
    n2 = n2/max(n2);
    % back to plotting
    figure(1)
    clf
    bar(x1, n1,'FaceColor','b','EdgeColor','none', 'FaceAlpha', .5);
    hold on
    bar(x2,n2, 'Facecolor', 'r', 'EdgeColor','none', 'FaceAlpha', .5);
    title(['Eigenmode Distribution ', num2str(filters(f))])
    xlabel('Eigenvalue')
    ylabel('Number of Modes')
    legend('Recovered Distribution', 'True Distribution')
    saveas(gca, ['/Users/tnl/Desktop/MATLAB/RAM/anesthesia/dci/filter_test/LPF_', num2str(filters(f)), '.jpg'], 'jpg');
end

