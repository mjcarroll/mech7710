%% Homework 4 - Problem 3
% Design a “Navigation” type Kalman filter to estimate the states 
% [East, North,Radar_Bias, Psi, Gyro_Bias]. (Note: this is a non-linear 
% problem that requires an Extended Kalman Filter (EKF) to do correctly. 
% However we can solve the problem in one of two ways: (i) linearize the 
% equations about the nominal operating point and produce a constant A 
% matrix for that operating point, or (ii) simply update the A matrix 
% at every time step with our measurements or estimates. Download the data
% hw4_3 from the website and run the filter sampled at 5 Hz.

clear; clc; close all;
load hw4_3.txt;

