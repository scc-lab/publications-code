% Rushikesh Kamalapurkar
% ADP Network



%% Initialization
test = 0;
if test==1
    clearvars -except Wc*HC*
    test = 1;
    tf=7;
else 
    clear all
    test = 0;
    tf=30;
end
close all
clc

global XX_iC; % Data store memory for storing the values of x
global U_iC; % Data store memory for storing the values of u
global SIGMATHETAS_iC; % Data store for storing the values of sigth,
           % made of sigth matrices stacked in a row
global xDmGUS_iC; % Data store memory for storing the values of xDot-G*u

global EC;
global MUC;
global UC;
fs = 1000; % sampling frequency
EC=cell(1,fs*tf+1);
MUC=cell(1,fs*tf+1);
UC=cell(1,fs*tf+1);


%% Problem Definition
RC={0.1;0.1;0.1;0.1;0.1};
QC={10;10;10;10;10};
n=1; % Number of States
m=[1;1;1;1;1]; % Number of controls
Nset=[1;2;3;4;5]; % Agent set


% Graph parameters
A = [0 1 0 0 0;
     1 0 0 0 0;
     0 0 0 0 0;
     0 0 1 0 0
     0 0 0 1 0];
A0 = [1 0 0 0 0;
      0 0 0 0 0;
      0 0 1 0 0;
      0 0 0 0 0
      0 0 0 0 0];

A=A(1:5,1:5);
A0=A0(1:5,1:5);
RC=RC(1:5);
QC=QC(1:5);
m=m(1:5);

Nset=Nset(1:5);

D = diag(sum(A,2));
calL = D-A;
[SC,ASC,A0SC,~,NmC]=ADPCLNNTSGenerateSubgraphs(A,A0,Nset);

N=length(Nset); % Number of agents
%%  Control Gains
L=[4;4;4;5;8]; % Number of basis functions

% ADP
etac1 = [0.1;0.1;0.1;0.1;0.1];
etac2 = [10;10;10;10;10];
etaa1 = [5;5;5;5;5];
etaa2 = [0.1;0.1;0.1;0.1;0.1];
beta = [0.3;0.3;0.3;0.3;0.3];
nu = [0.005;0.005;0.005;0.005;0.005];
GammaBar = [10000;10000;10000;10000;10000];

xDi0C={0.75;0;1;0;0};
xDijC=  {  0 0.5    0     0  0;
        -0.5   0    0     0  0;
           0   0    0     0  0;
           0   0 -0.5     0  0;
           0   0    0  -0.5  0};
xDi0CPLOT={0.75;0.25;1;0.5;0};
etac1=etac1(1:5);
etac2=etac2(1:5);
etaa1=etaa1(1:5);
etaa2=etaa2(1:5);
beta=beta(1:5);
nu=nu(1:5);
GammaBar=GammaBar(1:5);
xDi0C=xDi0C(1:5);
xDijC=xDijC(1:5,1:5);
xDi0CPLOT=xDi0CPLOT(1:5);
%% System ID
model=0;
P=[2;2;2;2;2];
k = [500;500;500;500;500];
GammaThetaC={1*eye(P(1));1*eye(P(2));1*eye(P(3));1*eye(P(4));1*eye(P(5))};
M_theta = [10;10;10;10;10]; % Number of recorded points
kTheta=[30;30;25;20;30];
fs=1000;
NN = 5;                 % Order of polynomial fit
F = 9;                % Window length
[b,g] = sgolay(NN,F);   % Calculate S-G coefficients
XX_iC=cell(N,1);
U_iC=cell(N,1);
SIGMATHETAS_iC=cell(N,1);
xDmGUS_iC=cell(N,1);
for i=1:N
    XX_i=zeros(n,F);
    U_i=zeros(m(i),(F+1)/2);
    SIGMATHETAS_i=zeros(P(i),M_theta(i));
    xDmGUS_i=zeros(n,M_theta(i));
    
    XX_iC{i}=XX_i;
    U_iC{i}=U_i;
    SIGMATHETAS_iC{i}=SIGMATHETAS_i;
    xDmGUS_iC{i}=xDmGUS_i;
end


k=k(1:5);
GammaThetaC=GammaThetaC(1:5);
M_theta=M_theta(1:5);
kTheta=kTheta(1:5);

%% Initial Conditions stack
xHC={0;0;0;0;0};
xC={2;2;2;2;2};
thetaHC={0*ones(P(1),n);0*ones(P(2),n);0*ones(P(3),n);0*ones(P(4),n);0*ones(P(5),n)};
vecthetaHC=cell(N,1);
for i=1:N
    vecthetaHC{i}=reshape(thetaHC{i},P(i)*n,1);
end
% thetaHC={[0;1];[0;1];0*ones(P(3),n);0*ones(P(4),n);0*ones(P(5),n)};
WcHC={1*ones(L(1),1);1*ones(L(2),1);1*ones(L(3),1);1*ones(L(4),1);3*ones(L(5),1)};
WaHC=WcHC;
GammaC={reshape(500*eye(L(1)),L(1)^2,1);
        reshape(500*eye(L(2)),L(2)^2,1);
        reshape(500*eye(L(3)),L(3)^2,1);
        reshape(500*eye(L(4)),L(4)^2,1);
        reshape(500*eye(L(5)),L(5)^2,1)};

if test==1
    etaa1=0;etaa2=0;etaa3=0;etaa4=0;etaa5=0;
    etac1=0;etac2=0;etac3=0;etac4=0;etac5=0;
    WcHC=WcHC_copy;
    WaHC=WaHC_copy;
end

L=L(1:5);
P=P(1:5);
xHC=xHC(1:5);
xC=xC(1:5);
thetaHC=thetaHC(1:5);
WcHC=WcHC(1:5);
WaHC=WaHC(1:5);
GammaC=GammaC(1:5);

% Initial vector
Zdim=[n*N;n*N;n*sum(P);sum(L);sum(L);sum(L.^2)];
Z0=[cell2mat(xC);cell2mat(xHC);...
    cell2mat(vecthetaHC);...
    cell2mat(WcHC);cell2mat(WaHC);...
    cell2mat(GammaC)];

%% BE extrapolation data points
M=zeros(N,1);
e_histC = cell(1,N);
x_histC = cell(1,N);
x0_histC = cell(1,N);
for i=1:N
    Si=SC{i};
    ei_hist=linspace(-1,1,5);
    ejxi_hist=linspace(-1,1,3);
    temp1=[kron(ones(1,n),ei_hist) kron(ones(1,length(Si)*n),ejxi_hist)];
    tempsize=[length(ei_hist)*ones(1,n) length(ejxi_hist)*ones(1,length(Si)*n)];
    temp2=mat2cell(temp1,1,tempsize);
    M(i)=length(ei_hist)^n*(n*length(ejxi_hist))^(length(Si));
    dimeSi=n*(length(Si)+1);
    [temp3{1:dimeSi}]=ndgrid(temp2{:});
    temp4=[];
    for j=1:dimeSi
        temp4=[temp4 temp3{j}(:)]; %#ok<AGROW>
    end
    E_i_histC=mat2cell(reshape(temp4,dimeSi*M(i),1),dimeSi*ones(M(i),1));
    ASi=ASC{i};
    A0Si=A0SC{i};
    xDijSiC=xDijC(Si,Si);
    xDi0SiC=xDi0C(Si);
    x_i_histC=cell(size(E_i_histC));
    e_i_histC=cell(size(E_i_histC));
    for j=1:M(i)
        ESik = E_i_histC{j};
        [xSik, x0k] = ADPCLNNTSGenStateFromError(n,ASi,A0Si,ESik,xDijSiC,xDi0SiC);
        eSik = ESik(1:end-n);
        x_i_histC{j,1}=xSik;
        e_i_histC{j,1}=eSik;
    end
    e_histC{i}=e_i_histC;
    x_histC{i}=x_i_histC;
    x0_histC{i}=x0k;
end
% t=1
% ADPCLNNTSDynamics(t,Z0,Zdim,A,A0,Nset,P,L,M,e_histC,x_histC,x0_histC,xDijC,xDi0C,n,m,RC,QC,etac1,etac2,etaa1,etaa2,beta,nu,GammaBar,k,GammaThetaC,kTheta,fs,F,g,M_theta,N);
%% Integration
% options = odeset('OutputFcn',@odeplot,'OutputSel',[1:10],'RelTol',1e-3,'AbsTol',1e-3);
% figure(3)
try
    Zout = ode1(@(t,Z)ADPCLNNTSDynamics(t,Z,Zdim,A,A0,Nset,P,L,M,e_histC,x_histC,x0_histC,xDijC,xDi0C,n,m,RC,QC,etac1,etac2,etaa1,etaa2,beta,nu,GammaBar,k,GammaThetaC,kTheta,fs,F,g,M_theta,N),[0:0.001:tf],Z0);
catch
    emailSimStatus('jane.doe@ufl.edu','Integration Failed','Simulation Update')
end
%  
%%
for i=1:size(Zout,1)
ZC=mat2cell(Zout(i,:)',Zdim,1);
xC=mat2cell(ZC{1},n*ones(length(Nset),1));
xHC=mat2cell(ZC{2},n*ones(length(Nset),1));
vecthetaHC=mat2cell(ZC{3},n*sum(P),1);
WcHC=mat2cell(ZC{4},L);
WaHC=mat2cell(ZC{5},L);
GammaC=mat2cell(ZC{6},L.^2);
X(:,i)= cell2mat(xC);
XH(:,i)= cell2mat(xHC);
ThetaH(:,i)= cell2mat(vecthetaHC);
WcH(:,i)= ZC{4};
WaH(:,i)= ZC{5};
end

E = cell2mat(EC);
U = cell2mat(UC);
MU = cell2mat(MUC);

LCost=10*E.^2+0.1*MU.^2;
LCost=LCost';
Lintegral=zeros(size(LCost));
for i=1:size(LCost,1)
    Lintegral(i,:) = 0.001*trapz(LCost(1:i,:),1);
end

T=0:0.001:tf;
x0 = exp(-0.1*T);
X0=zeros(N,length(T));
for i=1:N
    X0(i,:)=xDi0CPLOT{i}+x0;
end
emailSimStatus('jane.doe@ufl.edu','Simulation Complete','Simulation Update');