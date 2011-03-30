%% Homework 4 - Problem 2
% Download the data hw4_2 from the website.  The data is in the form [t y].
%  Supose we want to design an estimator to estimate the bias in the
%  measurement y.  We believe that the bias (x) is constant.

clear; clc; close all;
load hw4_2.txt;
t = hw4_2(:,1); y = hw4_2(:,2);
plot(t,y);
sigma_v = 1;

%% Part A
% Run the Kalman filter estimator with Q_d = 0.  What happens at t> 100
% seconds.  Why?  Calculate the steady state Kalman gain L_ss.  Plot L(k).
% This is known as the filter "going to sleep."

R = sigma_v^2;
A = 1; C = 1; 
Q = 0;

%Preallocations
P_plus = zeros(size(t)); x_hat_minus = zeros(size(t)); L = zeros(size(t));
P_minus = zeros(size(t)); x_hat_plus = zeros(size(t)); E = zeros(size(t));

for ii = 1:length(hw4_2),
    % Compute Kalman Gain
    P_plus(ii) = P_minus(ii) * (P_minus(ii) + R)^-1;
    L(ii) = P_plus(ii) * inv(R);
    
    % Update estimator with y(k)
    E(ii) = y(ii) - x_hat_minus(ii);
    x_hat_plus(ii) = x_hat_minus(ii) + L(ii) * E(ii);
    
    % Update covariance
    P_plus(ii) = (1 - L(ii)) * P_minus(ii);
    
    % Propagate Forward
    x_hat_minus(ii+1) = A * x_hat_plus(ii);
    P_minus(ii+1) = A * P_plus(ii) * A' + Q;
end

plot(t,y,t,x_hat_plus,'g');
figure();
plot(t,L);

%% Part B

%% Part C
% Now filter the measurement using the first order LPF.

Q = [0.00001, 0.0001, 0.001, 0.01];
color = ['g','r','b', 'y'];

figure(); hold on;

for ii = 1:length(Q),
    numd = sqrt(Q(ii));
    dend = [1,-(1-sqrt(Q(ii)))];
    yf = filter(numd,dend,y,y(1));
    plot(t,yf,color(ii))
end