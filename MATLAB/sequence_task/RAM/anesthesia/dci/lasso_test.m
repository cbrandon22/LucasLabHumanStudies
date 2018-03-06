%% LASSO Test/Simulation
%
% test LASSO for autoregressive precesses
%
% @author JStiso Mrch 2017

%% Starting with not autoregressive
% my understanding, in this example, x has 100 obs of a 5D vector. Z is 100
% observations as well, and A is 5 coeff
x = randn(100,5);
A = [1;2;1;-3;1]; % Only two nonzero coefficients
Z = x*A + randn(100,1)*.1; % Small added noise

[B, fit] = lasso(x, Z);

lassoPlot(B, fit)
%% Simulate an AR process
clear

% 2 observtions, 1 dimenson, A = .7, .25
model = arima('Constant',0.5,'AR',{0.7,0.25},'Variance',.1);
Y = simulate(model, 50, 'numpaths', 5); % get 50 time points

%% Plug in to lasso

X = Y(1:end-1,:);
X_hat = Y(end,:);

[A, fitAR] = lasso(X', X_hat, 'alpha', 1);

lassoPlot(A, fitAR)

%% Multidimensional

p = 1; % order, leave s 1 for now, number of observations per electrode
m = 4; % size of x, number of electrodes in actuality
tp = 2000; % number of time points to simulate

A = rand(m,m*p);
err = randn(tp,m)*.1;
y = zeros(tp,m);
% add noise and calculate x_hat
for i = 2:tp
    y(i,:) = y(i-1,:)*A + err(i,:);
end
y = y(1:370,:);

X = y(1:end-1,:);
X_hat = y(end,:);

[w,Aest,Cest,SBC,FPE,th] = arfit(y,1,1);

[A, fitAR] = lasso(X', X_hat, 'alpha', 1);