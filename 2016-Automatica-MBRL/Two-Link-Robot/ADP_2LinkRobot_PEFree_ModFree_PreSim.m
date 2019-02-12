% Rushikesh Kamalapurkar
% ADP 2 Link PE Free Simulink File

%% Initialization
% clear all
close all
clc

%% Control Gains
% Problem Definition
tf = 50;
R=1;
Q=[10 0 0 0; 0 10 0 0; 0 0 1 0; 0 0 0 1];
pecutoff = 0.4;
L=10;

typ = 1; % 1 = ADP OR 2 = Compare with optimal with constant weights.
model = 0; % 1 = Model-based OR 2 = Model-free
projection = 1; % 1 = projection-based update law for WaH OR 0.
learning = 1; % 1 = Concurrent Learning APD OR 2 = ADP with PE.

% % Gains for model = 1, projection = 0.
% etac1 = 0.5;
% etac2 = 64;
% etaa1 = 29.5;
% etaa2 = 0.01;
% beta = 0.01;

etac1 = 0.5*(projection==0) + 1*(projection==1);
etac2 = 10*(projection==0) + 20*(projection==1);
etaa1 = 5*(projection==0) + 0.1*(projection==1);
etaa2 = 0.001*(projection==0);
beta = 0.1;
v = 0.005;
GammaBar = 100000;

%% Initial Conditions stack
x0 = [1;1;0;0];
if typ == 1
    WcH0 = [5 5 0 0 0 0 25 0 2 2]';
else
      WcH0 = [26.5474045002859,0.640217153410971,3.53028402910059,3.46408279633458,12.3465599032095,1.16819734143584,44.6741606476645,4.06301554513049,4.67054605755498,0.149036117452653];
      WcH0 = [24.6297684659567,1.18840698913485,2.23807515028947,2.63543366699068,1.16589636116696,0.925620360472737,44.2298481912473,11.2724861339746,3.80631458802382,0.112951648575195];
end
WaH0 = 1*WcH0;
Gamma0 = 100*eye(L);

%% System ID
xH0 = 0*x0;
theta0 = 1*[1; 1; 1; 1];

z0 = [x0;xH0;theta0];

k = 10;
GammaTheta = diag([95 20 170 20])*(projection==0)+...
    diag([90 40 160 40])*(projection==1);
kTheta=1.1*(projection==0) + 1.1*(projection==1);
M = 20;
fs=1000;

NN = 5;                 % Order of polynomial fit
F = 9;                % Window length
[b,g] = sgolay(NN,F);   % Calculate S-G coefficients

%% Concurrent learning - Create database
Xc = linspace(-0.5,0.5,3);
N = length(Xc)^4;
index=1;
XQX=zeros(N,1);
PHIPF1=zeros(L,N);
PHIPY=[];
GG=[];
p1=3.473;
p2=.196;
p3=.242;
for i=1:length(Xc)
    for ii=1:length(Xc)
        for iii=1:length(Xc)
            for iiii=1:length(Xc)
                xc1=Xc(i);xc2=Xc(ii);xDc1=Xc(iii);xDc2=Xc(iiii);
                xc=[xc1 xc2 xDc1 xDc2]';
                Minvc = (1/(p2^2 - p1*p2 + p3^2*(cos(xc2))^2))*...
                        [-p2               p2+p3*cos(xc2)  ;...
                          p2+p3*cos(xc2) -p1-2*p3*cos(xc2)];
                Vmc = [ -p3*sin(xc2)*xDc2  , -p3*sin(xc2)*(xDc1+xDc2);
                        p3*sin(xc2)*xDc1  , 0                 ];
                Gc=[0 0; 0 0; Minvc]; 
                f1c=[xDc1; xDc2;-1*Minvc*Vmc*[xDc1; xDc2]];
                f2c=[Minvc Minvc]*diag([xDc1; xDc2;tanh(xDc1);tanh(xDc2)]);
                Yc=[0 0 0 0;0 0 0 0;-f2c];
%                 phic=...
%                       [xc(1)*xc(3);...
%                        xc(2)*xc(4);...
%                        xc(3)*xc(2);...
%                        xc(4)*xc(1);...
%                        xc(1)*xc(2);...
%                        xc(4)*xc(3);...
%                        xc(1)^2;...
%                        xc(2)^2;...
%                        xc(3)^2;...
%                        xc(4)^2];
                phi_pc=...
                      [xc(3)   0       xc(1)   0;...
                       0       xc(4)   0       xc(2);...
                       0       xc(3)   xc(2)   0;... 
                       xc(4)   0       0       xc(1);...
                       xc(2)   xc(1)   0       0;...
                       0       0       xc(4)   xc(3);...
                       2*xc(1) 0       0       0;...   
                       0       2*xc(2) 0       0;...
                       0       0       2*xc(3) 0;...
                       0       0       0       2*xc(4)];...
                XQX(index)=xc'*Q*xc;
                PHIPF1(:,index)=phi_pc*f1c;
                PHIPY=[PHIPY;phi_pc*Yc];
                GG=[GG;phi_pc*Gc*(R\Gc')*phi_pc'];
                index = index+1;
            end
        end
    end
end
