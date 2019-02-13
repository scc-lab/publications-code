% Rushikesh Kamalapurkar
% ADP 2 Link tracking PE Free Simulink File

%% Initialization
clear all
close all
clc

%% Control Gains
% Problem Definition
tf = 10;
R=1;
Q=[10 0 ; 0 10];
n=2; % number of states
m=1; % number of controls
model = 0; % 1 = Model-based OR 0 = Model-free

% Gains for model = 1, projection = 0.
s=0.7;
etac1 = 0.001;
etac2 = 2;
etaa1 = 2;
etaa2 = 0.0010;
beta = 0.01;
v = 1;
GammaBar = 1000;

% s=0.7;
% etac1 = 0.1;
% etac2 = 1.0;
% etaa1 = 3.0;
% etaa2 = 0.0001;
% beta = 0.01;
% v = 1;
% GammaBar = 10000;

MovePoints = 1;
ar = 3;
%% Initial Conditions stack
[~,~,L]=ADPSEStaF1LCLTr_Basis(zeros(2*n,1),s,n);
x0 = [0; 0];
WcH0 = 0.025*ones(L,1);
WaH0 = 0.025*ones(L,1);
Gamma0 = 50*eye(L);

freq = unifrnd(80,100,4,100);
Phase = unifrnd(0,2*pi,4,100);

% T1 = 0.1:pi/30:100;
% T2 = 0.1:exp(1)/25:100;
% T3 = 0.1:sqrt(5)/25:100;
% T4 = 0.1:sqrt(7)/30:100;

% freq = [T1(1:200);T1(1:200);T1(1:200);T1(1:200)];
% Phase = zeros(4,200);
% excite = [];
% for t=0:0.001:20
%     excite = [excite ADPSEStaF1LCLTr_ExtTraj(t,freq,Phase)];
% end
% XC = unifrnd(-ar*s,ar*s,2*n,1000);
%% System ID
xH0 = x0;
k = 500;
[~,outlayer]=ADPSEStaF1LCLTr_SIDBasis(x0);
GammaTheta = 1*eye(outlayer);
theta0=0*ones(outlayer,n);
M = 10; % Number of recorded points
kTheta=M*2;
fs=1000;
NN = 5;                 % Order of polynomial fit
F = 9;                % Window length
[b,g] = sgolay(NN,F);   % Calculate S-G coefficients