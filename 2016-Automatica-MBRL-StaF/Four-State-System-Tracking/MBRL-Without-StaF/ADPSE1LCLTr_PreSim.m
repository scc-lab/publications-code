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
typ = 1; % 1 = ADP OR 2 = Compare with optimal with constant weights.
model = 0; % 1 = Model-based OR 0 = Model-free
projection = 0; % 1 = projection-based update law for WaH OR 0.
basis = 1;

% Gains for learning=1, model = 1, projection = 0.
% etac1 = 0.1;
% etac2 = 1;
% etaa1 = 0.5;
% etaa2 = 0.01;
% beta =0.03;
% v = 0.005;
% GammaBar = 10000;
etac1 = 0.1;
etac2 = 2.5;
etaa1 = 1;
etaa2 = 0.01;
beta = 0.3;
v = 0.005;
GammaBar = 100000;

%% Initial Conditions stack
[~,L]=ADPSE1LCLTr_Basis(zeros(2*n,1),basis);
x0 = [0; 0];
if typ == 1
WcH0 = 5*ones(L,1);
WaH0 = 3*ones(L,1);
else
% WcH0 = [10.0030    4.4094    0.9064   -0.5879    1.0977   -0.0859   -0.0064    0.0973   -0.2685]';

WcH0 = [5    2.5   0    0   0   0    0.2]';
WaH0 = WcH0;
end
Gamma0 = 5000*eye(L);

%% System ID
xH0 = 0*x0;
p=25; % number of neurons
bias=0; % 0 for NN without bias, 1 for nn with bias.
k = 500;
Sidbasis=15;
if model==1
    Sidbasis=15;
end
inlayer=(n+1)*(bias==1)+n*(bias==0);
Ynn = -2.5*(0.5)+2.5*rand(inlayer,p);
[~,outlayer]=ADPSE1LCLTr_SIDBasis(x0,Sidbasis,bias,Ynn);
GammaTheta = 1*eye(outlayer);
if Sidbasis==1
    theta0=-8+16*rand(outlayer,n);
else
    theta0=0*rand(outlayer,n);
end
M = 10; % Number of recorded points
kTheta=M*2;
fs=1000;
NN = 5;                 % Order of polynomial fit
F = 9;                % Window length
[b,g] = sgolay(NN,F);   % Calculate S-G coefficients

%% Concurrent learning - Create database
Ec = linspace(-2,2,3);
tc = linspace(0.1,2*pi,12);
% N = length(tc)*length(Ec)^n;
N = length(Ec)^(2*n);
index=1;
QQ=zeros(N,1);
SIGPF1=zeros(L,N);
SIGP=zeros(N*L,2*n);
SIGPGD=zeros(N*L,n);
SIGTH=zeros(N*outlayer,N);
SIGTHD=zeros(N*outlayer,N);
GSIGMA=zeros(N*L,L);
for i=1:length(Ec)
    for ii=1:length(Ec)
        for iii=1:length(Ec)
%         for iii=1:length(tc)
            for iiii=1:length(Ec)
            xdc = [Ec(iii);Ec(iiii)];
%             xdc = [2*sin(tc(iii)); 2*cos(tc(iii))+2*sin(tc(iii))];
            hdc = [xdc(2)-xdc(1);xdc(2)-2*xdc(1)];
            ec1=Ec(i);ec2=Ec(ii);
            ec=[ec1 ec2]';
            xc=ec+xdc;
            Zc=[ec;xdc];
            gc=[0; (cos(2*xc(1))+2)];
            gdc = [0; (cos(2*xdc(1))+2)]; 
            gplusdc = [0 1/(cos(2*xdc(1))+2)];
            Gc = [gc;zeros(size(gc))];
            F1c = [-hdc+gc*gplusdc*hdc;hdc];
            [sigthc,~]=ADPSE1LCLTr_SIDBasis(xc,Sidbasis,bias,Ynn);
            [sigthdc,~]=ADPSE1LCLTr_SIDBasis(xdc,Sidbasis,bias,Ynn);
            [sig_pc,~]=ADPSE1LCLTr_Basis(Zc,basis);
            QQ(index)=ec'*Q*ec;
            SIGPF1(:,index)=sig_pc*F1c;
            SIGP(index*L-(L-1):index*L,:)=sig_pc;
            SIGPGD(index*L-(L-1):index*L,:)=sig_pc*[gc*gplusdc;zeros(n)];
            GSIGMA(index*L-(L-1):index*L,:)=sig_pc*Gc*(R\Gc')*sig_pc';
            SIGTH(index*outlayer-(outlayer-1):index*outlayer,index)=sigthc;
            SIGTHD(index*outlayer-(outlayer-1):index*outlayer,index)=sigthdc;
            index = index+1;
            end
        end
    end
end
