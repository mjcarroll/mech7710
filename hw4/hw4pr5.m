%% Homework 4 - Problem 5 (Optional)
% Estimator of a non-minimum phase system.

%% Part A
% Find the transfer function G_yw(s) from w(s) to y(s) and identify the
% poles and zeros.

A = [0,0;-1,-2];
B = [1;1];
C = [0,1];
D = 0;

[num,den] = ss2tf(A,B,C,D);
T = tf(num,den);

[z,p,k] = zpkdata(T);
zeros = z{1}
poles = p{1}

tzero(A,B,C,D)

%% Part B
% Sketch a root-locus.
%
% My root locus shows that the pole will travel towards s=-1 zero, meaning
% that the gain will never cause the estimation error to attenuate any
% faster than e^-t.

zeros = [zeros,-zeros]
poles = [poles(:);-poles(:)]

mirror = zpk(zeros,poles,1)

rlocus(mirror)

%% Part C
% Use the solution of the steady-state Riccati equation
%
% No matter how small Q (or the ratio w/v), the P matrix will never go
% completely to zero.
q = 1./logspace(1,10,100);
for ii = 1:length(q)
    P = care([0,0;-1,-2],[1;1],q(ii)*eye(2),1);
    test(ii) = P(1,1);
end

semilogx(q,test)
%% Part D

sys = ss(A,B,C,D);

%% Part E
% The results basically show that any system with non-minimum phase will
% not do as well in the Kalman filter/estimator setup.  This would be any
% sort of system with a considerable lag in it, or an initially negative
% step response.  It is quite possible to encounter these systems in the
% future, so the limitations of the Kalman filter are important to see.
