% Rushikesh Kamalapurkar
% ADP 2 Link tracking PE Free Simulink File

%% Initialization
clear all
close all
clc

%% Control Gains
% Problem Definition
tf = 10;
R1=1;
Q1=[10 0 0; 0 10 0; 0 0 10] ;
R2=1;
Q2=[10 0 0; 0 10 0; 0 0 10] ;
R3=1;
Q3=[10 0 0; 0 10 0; 0 0 10] ;
R4=1;
Q4=[10 0 0; 0 10 0; 0 0 10] ;
R5=1;
Q5=[10 0 0; 0 10 0; 0 0 10] ;
n=3; % number of states
m=2; % number of controls
model = 1; % 1 = Model-based OR 0 = Model-free

xd31=[-0.5;-0.5;0];
xd21=[0.5;-0.5;0];
xd23=xd21-xd31;
xd40=[0.5;0.5;0];
xd45=[1;0;0];
xd54=[-1;0;0];
% Gains for model = 1, projection = 0.
s1=0.25;
etac11 = 0.05;
etac21 = 5;
etaa11 = 1;
etaa21 = 2;
beta1 = 0.3;
v1 = 1;
GammaBar = 1000;

s3=0.25;
etac13 = 0.001;
etac23 = 5;
etaa13 = 1;
etaa23 = 2;
beta3 = 0.3;
v3 = 1;
GammaBar = 1000;

s2=0.25;
etac12 = 0.001;
etac22 = 5;
etaa12 = 0.5;
etaa22 = 1.5;
beta2 = 0.3;
v2 = 1;
GammaBar = 1000;

s4=0.25;
etac14 = 0.001;
etac24 = 5;
etaa14 = 0.5;
etaa24 = 1;
beta4 = 0.2;
v4 = 1;
GammaBar = 1000;

s5=0.25;
etac15 = 0.001;
etac25 = 5;
etaa15 = 1;
etaa25 = 1;
beta5 = 0.1;
v5 = 1;
GammaBar = 1000;
MovePoints = 1;
%% Initial Conditions stack
[~,~,~,L1]=ADPSEMREMKNN_Basis1(zeros(2*n,1));
x01 = [0; 0; 0];
WcH01 = 5*ones(L1,1);
WaH01 = 2*ones(L1,1);
Gamma01 = 20*eye(L1);

[~,~,~,~,L3]=ADPSEMREMKNN_Basis3(zeros(3*n,1));
x03 = [0; 0; 0];
WcH03 = 2*ones(L3,1);
WaH03 = 2*ones(L3,1);
Gamma03 = 30*eye(L3);

[~,~,~,~,~,L2]=ADPSEMREMKNN_Basis2(zeros(4*n,1));
x02 = [0; 0; 0];
WcH02 = 2*ones(L2,1);
WaH02 = 2*ones(L2,1);
Gamma02 = 30*eye(L2);

[~,~,~,~,L4]=ADPSEMREMKNN_Basis4(zeros(3*n,1));
x04 = 0*[0.5; -0.5; 0];
WcH04 = 5*ones(L4,1);
WaH04 = 1*ones(L4,1);
Gamma04 = 50*eye(L4);

[~,~,~,~,L5]=ADPSEMREMKNN_Basis5(zeros(3*n,1));
x05 = 0*[-0.5; -0.5; 0];
WcH05 = 2*ones(L5,1);
WaH05 = 1*ones(L5,1);
Gamma05 = 100*eye(L5);