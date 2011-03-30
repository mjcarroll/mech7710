%% Homework 4 - Part 4
% Estimator for vehicle dynamics
clear; clc;

Fs = 10;  Ts = 1/Fs;

%% Part A
% Assuming that we can only measure the yaw rate, design a Kalman filter to
% do full state estimation.  Provide a unit step steer input and estimate
% both states.  On one page plot the actual states and estimated states.
% Where are the steady state poles of the estimator?

A = [-2.62, 12; -0.96, -2]; B = [14;1]; C = [1 0]; D = 0;
[Ad,Bd,Cd,Dd] = c2dm(A,B,C,D,Ts,'zoh');

plant_d = ss(Ad,[Bd Bd],Cd,Dd,Ts,...
    'inputname',{'u','w'}, ...
    'statename',{'yawangle','slipangle'});

% Use kalman to design the Kalman filter
Q = 0.001; R = 0.1^2;
[kalmf,L,P,M] = kalman(plant_d,Q,R);

% Construct the simulated system
% Add some augmented outputs so that I can get at the states
a_sim = Ad;
b_sim = [Bd, Bd, 0*Bd];
c_sim = [Cd; Cd; [0,1]];
d_sim = [0,0,1;0,0,0;0,0,0];

sys_sim = ss(a_sim,b_sim,c_sim,d_sim,Ts,'inputname',{'u','w','v'},...
    'outputname',{'y1','x1','x2'},...
    'statename',{'x1','x2'});
SimSys = parallel(sys_sim,kalmf,1,1,[],[]);
SimSys = feedback(SimSys,1,4,1,1);
SimSys = SimSys([1 2 3 4 5 6],[1 2 3]);

Q = 0.001; R = 0.1^2;
t = [0:0.1:20]';
u = [zeros(1,100),ones(1,length(t)-100)]';
n = length(t);
randn('seed',0);
w = sqrt(Q) * randn(n,1);
v = sqrt(R) * randn(n,1);

out = lsim(SimSys,[w,v,u]);

figure();
plot(t,out(:,1),t,out(:,4),'g');
figure();
subplot(2,2,1); plot(t,out(:,2)); title('x1')
subplot(2,2,2); plot(t,out(:,5)); title('x1_e')
subplot(2,2,3); plot(t,out(:,3)); title('x2')
subplot(2,2,4); plot(t,out(:,6)); title('x2_e')

%% Part B
% Part A now with a changed center of gravity
A = [-2.62, 12; -0.96, -2]; B = [14;1]; C = [1 0]; D = 0;
[Ad,Bd,Cd,Dd] = c2dm(A,B,C,D,Ts,'zoh');

plant_d = ss(Ad,[Bd Bd],Cd,Dd,Ts,...
    'inputname',{'u','w'}, ...
    'statename',{'yawangle','slipangle'});

% Use kalman to design the Kalman filter
Q = 0.001; R = 0.1^2;
[kalmf,L,P,M] = kalman(plant_d,Q,R);

A = [-2.42 4; -0.99, -2]; B = [18;1]; C = [1 0]; D = 0;
[Ad,Bd,Cd,Dd] = c2dm(A,B,C,D,Ts,'zoh');

a_sim = Ad;
b_sim = [Bd, Bd, 0*Bd];
c_sim = [Cd; Cd; [0,1]];
d_sim = [0,0,1;0,0,0;0,0,0];

sys_sim = ss(a_sim,b_sim,c_sim,d_sim,Ts,'inputname',{'u','w','v'},...
    'outputname',{'y1','x1','x2'},...
    'statename',{'x1','x2'});
SimSys = parallel(sys_sim,kalmf,1,1,[],[]);
SimSys = feedback(SimSys,1,4,1,1);
SimSys = SimSys([1 2 3 4 5 6],[1 2 3]);

t = [0:0.1:20]';
u = [zeros(1,100),ones(1,length(t)-100)]';
n = length(t);
randn('seed',0);
w = sqrt(Q) * randn(n,1);
v = sqrt(R) * randn(n,1);

out = lsim(SimSys,[w,v,u]);

figure(3);
plot(t,out(:,1),t,out(:,4),'g');
figure(4);
subplot(2,2,1); plot(t,out(:,2)); title('x1')
subplot(2,2,2); plot(t,out(:,5)); title('x1_e')
subplot(2,2,3); plot(t,out(:,3)); title('x2')
subplot(2,2,4); plot(t,out(:,6)); title('x2_e')

%% Part C
% The new R will be <LATEX HERE>

%% Part D
% Redo part A with the noisy slip angle

A = [-2.62, 12; -0.96, -2]; B = [14;1]; C = [1 0; 0 1]; D = [0;0];
[Ad,Bd,Cd,Dd] = c2dm(A,B,C,D,Ts,'zoh');

plant_d = ss(Ad,[Bd Bd],Cd,0,Ts,...
    'inputname',{'u','w'}, ...
    'statename',{'x1','x2'});

% Use kalman to design the Kalman filter
Q = 1; R = [0.1^2,0; 0 0.5^2];
[kalmf,L,P,M] = kalman(plant_d,Q,R,0);

% Construct the simulated system
% Add some augmented outputs so that I can get at the states
a_sim = Ad;
b_sim = [Bd, Bd, 0*Bd,0*Bd];
c_sim = [Cd;Cd];
d_sim = [0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];

sys_sim = ss(a_sim,b_sim,c_sim,d_sim,Ts,...
    'inputname',{'u','w','v1','v2'},...
    'outputname',{'y1','y2','x1','x2'},...
    'statename',{'x1','x2'});
SimSys = parallel(sys_sim,kalmf,1,1,[],[]);
SimSys = feedback(SimSys,1,5,1,1);
SimSys = feedback(SimSys,1,6,2,1);
SimSys = SimSys([1 2 3 4 5 6 7 8],[1 2 3 4]);

t = [0:0.1:20]';
u = [zeros(1,100),ones(1,length(t)-100)]';
n = length(t);
randn('seed',0);
w = sqrt(Q) * randn(n,1);
v1 = sqrt(R(1,1)) * randn(n,1);
v2 = sqrt(R(2,2)) * randn(n,1);

out = lsim(SimSys,[w,v1,v2,u]);
figure(5);
subplot(2,1,1); plot(t,out(:,1),t,out(:,5),'g');
subplot(2,1,2); plot(t,out(:,2),t,out(:,6),'g');
figure(6);
subplot(2,2,1); plot(t,out(:,3)); title('x1')
subplot(2,2,2); plot(t,out(:,7)); title('x1_e')
subplot(2,2,3); plot(t,out(:,4)); title('x2')
subplot(2,2,4); plot(t,out(:,8)); title('x2_e')

mean(out(:,3) - out(:,7))
mean(out(:,4) - out(:,8))