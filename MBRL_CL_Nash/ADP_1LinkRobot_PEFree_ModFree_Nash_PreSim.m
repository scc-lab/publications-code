% Rushikesh Kamalapurkar
% ADP 1 Link PE Free Simulink File

%% Initialization
clear all
close all
clc

%% Control Gains
% Problem Definition
tf = 10;
n=2; % State Dimension
m=1; % Control Dimension
pecutoff = .4;
L=3;
model = 1; % 1 = Model-based OR 0 = Model-free
learning = 1; % 1 = Concurrent Learning APD OR 2 = ADP with PE.
test = 0;
% player 1 gains
Q_1=2*eye(2);
R_11=2;
R_21=1;
etac1_1 = 1; %1
etac2_1 = 1.5; %1
etaa1_1 = 10; %10
etaa2_1 = 0.1; %0.1
beta_1 = 3; %3
v_1 = 0.005; %0.005
GammaBar_1 = 10000;
% player 2 gains
Q_2=eye(2);
R_22=1;
R_12=2;
etac1_2 = 1; %1
etac2_2 = 1; %1
etaa1_2 = 10; %10
etaa2_2 = 0.1; %0.1
beta_2 = 3; %3
v_2 = 0.005; %0.005
GammaBar_2 = 10000;

%% Initial Conditions stack
x0 = [1;1];
% player 1
WcH0_1 = 3*[1 1 1]';
WaH0_1 = 1*WcH0_1;
Gamma0_1 = 100*eye(L);
% player 2
WcH0_2 = 3*[1 1 1]';
WaH0_2 = 1*WcH0_2;
Gamma0_2 = 100*eye(L);

%% System ID
xH0 = 0*x0;
theta0 = 0.5*ones(1,6)'; % actual theta: [1 -2 -.5 -1 .25 .25]
z0 = [x0;xH0;theta0];
k = 5;
GammaTheta = 20*diag([1 1 5 5 3 3]); % [1 1 2 5 3 3]
kTheta=1.5;
M = 30;
fs = 1000;

NN = 5;                 % Order of polynomial fit
F = 9;                  % Window length
[b,g] = sgolay(NN,F);   % Calculate S-G coefficients

%% Concurrent learning - Create database
% used for both players

% player 1
DY1 = 2; DY2 = 6; % Y is a 2 x 5 matrix 
Xc = linspace(-2,2,5);
N = length(Xc)^n;
index=1;
XQX1=zeros(N,1);
XQX2=zeros(N,1);
PHIPY=[];
GG1=[];
GG2=[];
GG12=[];
GG21=[];
for i=1:length(Xc)
    for ii=1:length(Xc)
        xc1=Xc(i);xc2=Xc(ii); xc = [xc1;xc2];
        Gc1 = [0 ; cos(2*xc1)+2];
        Gc2 = [0 ; sin(4*xc1^2)+2]; 
        Yc = [xc(2) xc(1) 0 0 0 0 ;0 0 xc1 xc2 xc2*(cos(2*xc1)+2)^2 xc2*(sin(4*xc1^2)+2)^2];
        phi_pc =[2*xc1 0; xc2 xc1; 0 2*xc2];
        XQX1(index)=xc'*Q_1*xc;
        XQX2(index)=xc'*Q_2*xc;
        PHIPY=[PHIPY;phi_pc*Yc];
        GG1=[GG1;phi_pc*Gc1*(R_11\Gc1')*phi_pc'];
        GG12=[GG12;phi_pc*Gc2*(1/R_22)*R_12*(1/R_22)*Gc2'*phi_pc'];
        GG21=[GG21;phi_pc*Gc1*(1/R_11)*R_21*(1/R_11)*Gc1'*phi_pc'];
        GG2=[GG2;phi_pc*Gc2*(R_22\Gc2')*phi_pc'];
        index = index+1;
    end
end

