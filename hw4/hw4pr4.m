%% Homework 4 - Part 4
% Estimator for vehicle dynamics
clear; clc;

Fs = 10;  Ts = 1/Fs;
A = [-2.62, 12; -0.96, -2]; B = [14;1]; C = [1 0]; D = 0;

plant = ss(A,B,C,0,...
    'inputname','u',...
    'outputname','y',...
    'statename',{'YawRate','SideSlip'});

%% Part A
% Assuming that we can only measure the yaw rate, design a Kalman filter to
% do full state estimation.  Provide a unit step steer input and estimate
% both states.  On one page plot the actual states and estimated states.
% Where are the steady state poles of the estimator?

[Ad,Bd,Cd,Dd] = c2dm(A,B,C,D,Ts,'zoh');

plant_d = ss(Ad,[Bd Bd],Cd,0,Ts,...
    'inputname',{'u','w'});

Q = 1; R = 0.1^2;

[kalmf,L,P,M] = kalman(plant_d,Q,R);
kalmf = kalmf(1,:);

% Construct the simulated system
a_sim = Ad;
b_sim = [Bd, Bd, 0*Bd];
c_sim = [Cd; Cd];
d_sim = [0,0,0;0,0,1];

sys_sim = ss(a_sim,b_sim,c_sim,d_sim,'inputname',{'u','w','v'},...
    'outputname',{'y','yv'});

sys = parallel(sys_sim,kalmf,1,1,[],[]);