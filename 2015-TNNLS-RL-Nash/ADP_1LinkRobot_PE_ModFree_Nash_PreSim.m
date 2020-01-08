% Rushikesh Kamalapurkar
% ADP 1 Link PE Free Simulink File

%% Initialization
clear all
close all
clc

%% Control Gains
% Problem Definition
tf = 20;
n=2; % State Dimension
m=1; % Control Dimension
pecutoff = .3;
L=3;
Lf=5;
model = 0; % 1 = Model-based OR 0 = Model-free
learning = 1; % 1 = Concurrent Learning APD OR 2 = ADP with PE.
test = 0;
%Derivative estimator gains
k=300;%500
alpha=200;%300
gamma=5;
beta1=0.2;
GammaDNN=0.1;
% % player 1 gains
% Q_1=2*eye(2);
% R_11=2;
% R_21=1;
% etac1_1 = 15; %15
% etaa1_1 = 10; %10
% etaa2_1 = 50; %50
% beta_1 = 0.03; %0.03
% v_1 = 0.001; %0.001
% % player 2 gains
% Q_2=eye(2);
% R_22=1;
% R_12=2;
% etac1_2 = 15; %15
% etaa1_2 = 10; %10
% etaa2_2 = 50; %50
% beta_2 = 0.03; %0.03
% v_2 = 0.001; %0.001

% player 1 gains
Q_1=2*eye(2);
R_11=2;
R_21=1;
etac1_1 = 25; %15
etaa1_1 = 0.1; %10
etaa2_1 = 20; %50
beta_1 = 0.03; %0.03
v_1 = 0.001; %0.001
% player 2 gains
Q_2=eye(2);
R_22=1;
R_12=2;
etac1_2 = 10; %15
etaa1_2 = 0.1; %10
etaa2_2 = 20; %50
beta_2 = 0.03; %0.03
v_2 = 0.001; %0.001

%% Initial Conditions stack
x0 = [3;-1];
xH0 = [0;0];
nu0 = [0;0]; 
WfH0 = -1+2*rand((Lf+1),n);
VfH0 = -1+2*rand(n,Lf);
% player 1
WcH0_1 = 3*[1 1 1]';
WaH0_1 = 1*WcH0_1;
Gamma0_1 = 5000*eye(L);
% player 2
WcH0_2 = 3*[1 1 1]';
WaH0_2 = 1*WcH0_2;
Gamma0_2 = 5000*eye(L);