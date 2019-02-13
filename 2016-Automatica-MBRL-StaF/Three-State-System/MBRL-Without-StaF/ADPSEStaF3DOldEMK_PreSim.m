% Rushikesh Kamalapurkar
% ADP 1 Link PE Free Simulink File

%% Initialization
clear all
close all
clc

%% Control Gains
% Problem Definition
tf = 5;
n=3; % State Dimension
m=2; % Control Dimension
R=eye(m);
Q=10*eye(n);
pecutoff = 0.6;

model = 0; % 1 = Model-based OR 0 = Model-free
learning = 1; % 1 = Concurrent Learning APD OR 2 = ADP with PE.
etac1 = 0.001;
etac2 = 20;
etaa1 = 17;
etaa2 = 0.001;
beta = 0.0001;
v = 0.5;
s=1;
dbar=0.0;
v2 = 2;
MovePoints = 0;
NormBas = 0;
NormDyn = 0;

L=9;
D=0;
% etac1 = 10;
% etac2 = 0;
% etaa1 = 10;
% etaa2 = 0;
% beta = 0;
% v = 0.005;
GammaBar = 100000;

%% Initial Conditions stack
x0 = [-0.5; 0.5; -0.250];
% WcH0 = 3*[1 1 1]';
WcH0=1*ones(L,1);
WaH0 = 0.85*WcH0;
Gamma0 = 100*eye(L);
% Gamma0 = 5000*eye(L);

% %% System ID
% xH0 = 0*x0;
% theta0 = 1*[1; 1; 1; 1];
% 
% z0 = [x0;xH0;theta0];
% % 
% k = 10;
% GammaTheta = diag([20 20 20 20]);
% kTheta=30;
% 
% % k = 10;
% % GammaTheta = diag([20 20 20 20]);
% % kTheta=1;
% 
% % k = 5;
% % GammaTheta = 10*diag([1 1 2 4 1 1]);
% % kTheta=1.2;
% 
% M = 30;
% fs = 1000;
% 
% NN = 5;                 % Order of polynomial fit
% F = 9;                % Window length
% [b,g] = sgolay(NN,F);   % Calculate S-G coefficients
% 
