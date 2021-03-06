% Rushikesh Kamalapurkar
% ADP 1 Link PE Free Simulink File

%% Initialization
clear all
close all
clc

%% Control Gains
% Problem Definition
tf = 5;
n=2; % State Dimension
m=1; % Control Dimension
R=1;
Q=[1 0 ; 0 1];
pecutoff = 0.6;
L=3;
model = 0; % 1 = Model-based OR 0 = Model-free
learning = 1; % 1 = Concurrent Learning APD OR 2 = ADP with PE.
etac1 = 0.001;
etac2 = 0.25;
etaa1 = 1.2;
etaa2 = 0.01;
beta = 0.003;
v = 0.05;
s=0.7;
dbar=0.00001;
v2 = 1;
MovePoints = 1;
NormBas = 0;
NormDyn = 0;
d1=s*[0;1];
d2=s*[0.8660;-0.5];
d3=s*[-0.8660;-0.5];
% d2=s*[0.6660;-0.5];
% d3=s*[-0.3660;-0.5];

% etac1 = 10;
% etac2 = 0;
% etaa1 = 10;
% etaa2 = 0;
% beta = 0;
% v = 0.005;
GammaBar = 100000;

%% Initial Conditions stack
x0 = [-1;1];
WcH0 = 0.4*[1 1 1]';
WaH0 = 0.7*WcH0;
Gamma0 = 250*eye(L);
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
%% Concurrent learning - Create database
DY1 = 2; DY2 = 4; % Y is a 2 x 4 matrix 
Xc = linspace(-2,2,3);
N = length(Xc)^n;
index=1;
XQX=zeros(N,1);
PHIPY=[];
GG=[];
for i=1:length(Xc)
    for ii=1:length(Xc)
        xc1=Xc(i);xc2=Xc(ii); xc = [xc1;xc2];
        Gc = [0 ; cos(2*xc1)+2]; 
        Yc = [xc1 xc2 0 0 ;0 0 xc1 xc2*(1-(cos(2*xc1)+2)^2)];
        phi_pc =[2*xc1 0; xc2 xc1; 0 2*xc2];
        XQX(index)=xc'*Q*xc;
        PHIPY=[PHIPY;phi_pc*Yc];
        GG=[GG;phi_pc*Gc*(R\Gc')*phi_pc'];
        index = index+1;
    end
end
