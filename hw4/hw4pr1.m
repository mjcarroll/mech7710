%% Homework 4 - Problem 1
% Kalman filter at its best - simulation.  Suppose we have a 2nd order
% system that we are regulating about zero (position and velocity) by
% wrapping an "optimal" control loop around the system.  The new dynamics
% of the continuous time system are given by the closed-loop A matrix:
%
% Suppose our measurement is simply position (C=[1 0]).  There is a white
% noise process disturbance (force, B_w = [0 1]^T) acting on the controlled
% system.

clear; clc; close all; format loose; format short;
Fs = 10; Ts = 1/Fs;
A_cl = [0, 1; -1, -1.4];
C = [1, 0];
Bw = [0, 1]';
D = [0,1];

%% Part A
% Simulate the controlled system with the disturbance force (1\sigma = 2)
% and a sampled sensor noise (1\sigma = 1) for 100 seconds at 10Hz.

sigma_q = 2;
sigma_r = 1;

t = 0:Ts:100;
w = sigma_q * randn(size(t));
v = sigma_r * randn(size(t));

[Ad,Bd,Cd,Dd] = c2dm(A_cl,Bw,C,D,Ts,'zoh');

dPlant = ss(Ad,[zeros(size(Bd)),Bd,zeros(size(Bd))],Cd,0,Ts,...
    'statename',{'position' 'velocity'},...
    'inputname',{'u','w','v'},...
    'outputname',{'y'});

[y_d,t_d,x_d] = lsim(dPlant,[zeros(size(t));w]);


%% Part B
% What is Q, Q_d and R_d?

Q = 2^2
R = sigma_r^2
% Bryson's Trick
S = [-A_cl, Bw*Q*Bw'; zeros(2), A_cl'];
C_bryson = expm(S*Ts);
Qd = C_bryson(3:4,3:4)' * C_bryson(1:2,3:4)
R_d = sigma_r^2/Ts

%% Part C
% Calculate the steady state Kalman gain for the system.  This can be done
% in one of many ways: iterate the Kalman filter until it converges,
% dlqe.m, dare.m, kalman.m, dlqr.m (+ predictor to current estimator
% trick), tec.  What is the steady state covariance of the estimates after
% the time update, P^{(-)}, as well as after the measurement update,
% P^{(+)}.  Where are the poles of the estimator?

%[M,P,Z,E] = lqed(A_cl, Bw, C, Q, R, Ts);
[M,P,Z,E] = dlqe(Ad, eye(2), Cd, Qd, R_d);

P_plus = P
P_minus = Ad * P_plus * Ad' + Qd
L = M
Poles = E
zplane([],E);

%% Part D
% Use the steady state Kalman filter to generate an estimate of the 2
% states over time.  Calculate the norm of the standard deviation of the
% errors for each state.

sys = ss(Ad,[zeros(size(Bw)),Bw],Cd,0,Ts,...
    'statename',{'position' 'velocity'},...
    'inputname',{'u','w'},...
    'outputname',{'y'});

[kalmf,L,P,M] = kalman(sys,Q,R);

