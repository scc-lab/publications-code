% Rushikesh Kamalapurkar
% ADP 2 Link tracking PE Free Simulink File

%% Initialization
clear all
close all
clc

%% Control Gains
% Problem Definition
tf = 20;
R=1;
Q=[10 0 ; 0 10];
n=2; % number of states
m=1; % number of controls
typ = 1; % 1 = ADP OR 2 = Compare with optimal with constant weights.
model = 0; % 1 = Model-based OR 0 = Model-free
projection = 1; % 1 = projection-based update law for WaH OR 0.
basis = 1;
learning=1;
pecutoff=0.5;

% Gains for learning=1, model = 1, projection = 0.
if learning==1
    etac1 = 0.5;
    etac2 = 10;
    etaa1 = 10;
    etaa2 = 0.001;
    beta = 0.1;
    v = 0.005;
    GammaBar = 100000;
else
    model=1;
    etac1 = 1;
    etac2 = 0;
    etaa1 = 0.1;
    etaa2 = 0.001;
    beta = 0.1;
    v = 0.005;
    GammaBar = 100000;
end
% etac1 = 0.01;
% etac2 = 15;
% etaa1 = 1;
% etaa2 = 0.01;
% beta = 0.3;
% v = 0.005;
% GammaBar = 10000;

%% Initial Conditions stack
[~,L]=ADPSELSCLTr_Basis(zeros(2*n,1),basis);
x0 = [1; 1];
if typ == 1
WcH0 = 10*ones(L,1);
WaH0 = 10*ones(L,1);
else
% WcH0 = [5.03281594074094,-0.000407478452205025,-0.00130553984605153,2.42849166101254,0.0538748951810672,-0.201714909511405,0.253851050892307,0.0226046803934718,1.87904944703855]';
% WaH0 = [5.03281594074094,-0.000407478452205025,-0.00130553984605153,2.42849166101254,0.0538748951810672,-0.201714909511405,0.253851050892307,0.0226046803934718,1.87904944703855]';
WcH0 = [0.00981102774694368,1.25468716257868,1.90522222702086,-0.221489942266716,3.96914457427608,-0.434624884727441,0.220704157419489]';
WaH0 = [0.00981102774694368,1.25468716257868,1.90522222702086,-0.221489942266716,3.96914457427608,-0.434624884727441,0.220704157419489]';
end

Gamma0 = 1000*eye(L);

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
[~,outlayer]=ADPSELSCLTr_SIDBasis(x0,Sidbasis,bias,Ynn);
GammaTheta = 1*eye(outlayer);
if Sidbasis==1
    theta0=-8+16*rand(outlayer,n);
else
    theta0=0*ones(outlayer,n);
end
M = 10; % Number of recorded points
kTheta=M*1;
fs=1000;
NN = 5;                 % Order of polynomial fit
F = 9;                % Window length
[b,g] = sgolay(NN,F);   % Calculate S-G coefficients

%% Concurrent learning - Create database
Ec = linspace(-2,2,5);
tc = linspace(0.1,2*pi,11);
% N = length(Ec)^(2*n);
N = length(tc)*length(Ec)^(2);
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
%             for iiii=1:length(Ec)
            xdc = [2*sin(tc(iii)); 2*cos(tc(iii))+2*sin(tc(iii))];
            hdc = [xdc(2)-xdc(1);xdc(2)-2*xdc(1)];
            ec1=Ec(i);ec2=Ec(ii);
%             xdc1=Ec(iii);xdc2=Ec(iiii);
            ec=[ec1 ec2]';
%             xdc=[xdc1 xdc2]';
%             hdc = [xdc(2)-xdc(1);xdc(2)-2*xdc(1)];
            xc=ec+xdc;
            Zc=[ec;xdc];
            gc=[0; 1];
            gdc = [0; 1]; 
            gplusdc = [0 1];
            Gc = [gc;zeros(size(gc))];
            A=[-1 1; -0.5 -0.5];
            Fc=[A*xc-gdc*gplusdc*A*xdc;[0;0]];
            F1c = [-hdc+gc*gplusdc*hdc;hdc];
%             if model==1
                F1cm=F1c+Fc;
%             end
            [sigthc,~]=ADPSELSCLTr_SIDBasis(xc,Sidbasis,bias,Ynn);
            [sigthdc,~]=ADPSELSCLTr_SIDBasis(xdc,Sidbasis,bias,Ynn);
            [sig_pc,~]=ADPSELSCLTr_Basis(Zc,basis);
            QQ(index)=ec'*Q*ec;
            SIGPF1(:,index)=sig_pc*F1c;
            SIGPF1m(:,index)=sig_pc*F1cm;
            SIGP(index*L-(L-1):index*L,:)=sig_pc;
            SIGPGD(index*L-(L-1):index*L,:)=sig_pc*[gc*gplusdc;zeros(n)];
            GSIGMA(index*L-(L-1):index*L,:)=sig_pc*Gc*(R\Gc')*sig_pc';
            SIGTH(index*outlayer-(outlayer-1):index*outlayer,index)=sigthc;
            SIGTHD(index*outlayer-(outlayer-1):index*outlayer,index)=sigthdc;
            index = index+1;
%             end
        end
    end
end
