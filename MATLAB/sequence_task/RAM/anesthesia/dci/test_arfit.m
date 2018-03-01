%% Simulate AR process in order to test AR fit
% Simulate and AR process, and check if arfit can give back original
% varibles
% X_hat_t = Ax_t-1 + e

%% Make random coefficient matrices and vectors

p = 1; % order, leave s 1 for now, number of observations per electrode
m = 1; % size of x, number of electrodes in actuality
tp = 2000; % number of time points to simulate
lag = 20; % dont use the fist few timepoints when checking

A = randn(m,p*m);
err = randn(tp,m);
x = zeros(tp,m);
% add noise and calculate x_hat
for i = 2:tp
    x(i,:) = x(i-1,:)*A + err(i,:);
end

[w,Aest,Cest,SBC,FPE,th] = arfit(x(:,:),1,1);

plot(x)

% get error between eigen vals
eig_s = eig(A);
eig_est = eig(Aest);
rmse = sqrt((eig_s-eig_est).^2)

% maybe x is nonstationary? or order 1 is not good?

%% Try higher orders

p = 2; % order, number of observations per electrode
m = 2; % size of x, number of electrodes in actuality
tp = 2000; % number of time points to simulate

A = rand(m,p*m);
err = randn(tp,m);
x = zeros(tp,m);
% add noise and calculate x_hat
for i = 2:tp
    for j = 1:p
        x(i,:) = x(i,:) + x(i-1,:)*A(:,1+((j-1)*m):j*m);
    end
    x(i,:) = err(i,:);
end

[w,Aest,Cest,SBC,FPE,th] = arfit(x,p,p);

plot(x)

% get error between eigen vals
eig_s = eig([A; eye((p-1)*m) zeros((p-1)*m)]);
eig_est = eig([Aest; eye((p-1)*m) zeros((p-1)*m)]);
rmse = sqrt((eig_s-eig_est).^2)

