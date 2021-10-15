%% Paper title "Safe Model-Based Reinforcement Learning for Systems with Parametric Uncertainties" 
%% arXiv:2007.12666
%% Coded by S M Nahid Mahmud, MS Grad Student, Oklahoma State University.
%% nahid.mahmud@okstate.edu
%% Four-state dynamical system (Section 7.2 of the journal paper)

clear all
clc
tic
%% PreSim

a1 = -7; A1 = 5; %Barrier Boundaries for x1, x1-coordinate values will be remained within (a1,A1) 
a2 = -7; A2 = 5; %Barrier Boundaries for x2, x2-coordinate values will be remained within (a2,A2) 
a3 = -5; A3 = 7; %Barrier Boundaries for x3, x3-coordinate values will be remained within (a3,A3) 
a4 = -5; A4 = 7; %Barrier Boundaries for x4, x4-coordinate values will be remained within (a4,A4) 

P.n = 4; %number of states
P.m = 1; %number of dimension
P.l = 4; %number of unknown parameters
P.x0 = [-5;-5;5;5]; % Initial x-coordinate values, it has to be within (a,A).

% known parameters, 
p1 = 3.473;
p2 = 0.196;
p3 = 0.242;

%transforming initial x-coordinates to initial s-coordinates. 
s10 = (log ((A1/a1)*((a1-P.x0(1))/(A1-P.x0(1)))));
s20 = (log ((A2/a2)*((a2-P.x0(2))/(A2-P.x0(2)))));
s30 = (log ((A3/a3)*((a3-P.x0(3))/(A3-P.x0(3)))));
s40 = (log ((A4/a4)*((a4-P.x0(4))/(A4-P.x0(4)))));

P.s0 = [s10;s20;s30;s40]; % Initial s-coordinate values
P.W = 10; % number of the Actor Weight vector elements as wellas Critic Weight vector elements.
P.G = 10;  % number of the Gamma vector elements. 

P.theta = [5.3;1.1;8.45;2.35]; % True Value of the unknown parameters 
P.theta0= [5;5;5;5]; % Picked initial value of the unknown parameters set by the user.
P.theta_tilde0 = P.theta-P.theta0; % Errors between the value of the unknown parameters and initial picked up value of the unknown parameters. 
P.IG = [0;0;0;0]; % Initial value of the integrated G function.  
P.G1 = 4;  % number of the elements in G vector.
P.G0 = eye(P.G);


P.Wa0 = [60;2;2;2;2;2;40;2;2;2]; % Initial value of the Actor weight vector.  
P.Wc0 = [60;2;2;2;2;2;40;2;2;2]; % Initial value of the Critic weight vector.  
P.Cost0=0; % Initial integrated cost. 

%%Initial values for s = 1-4, thetahat= 5-8, scriptY = 9-24, scriptYf =
%25-40, scriptXf = 41-44, WaH = 45:54, WcH = 55:64, Gamma =  65:164, 
% x = 165:168, Int_G = 169:172, Cost = 173.

% ThetaHat denotes Estimated value of the unknown parameters, WaH denotes Estimated Actor weight, WcH denotes Estimated Critic weight,
% Int_G denotes the integration of G function. 
P.z0 = [P.s0;P.theta0;zeros(P.n*P.l+P.l^2+P.l,1);P.Wa0;P.Wc0;P.G0(:);P.x0;P.IG;P.Cost0];


%% Integration
%ODE45 function
options = odeset('OutputFcn',@odeplot,'OutputSel',1:P.n);
[t,z] = ode45(@(t,z) closedLoopDynamics(t,z,P),[0 60],P.z0,options) ;%ODE45 function

%% Plot

% State trajectories of the s-coordinates
figure(1)
plot(t,z(:,1:P.n),'Linewidth',2);
xlabel('Time, t(s)','interpreter','latex')
ylabel('Transformed States, s (t)','interpreter','latex')
set(gca,'FontSize',10)
l=legend('$s_{1}$(t)','$s_{2}$(t)','$s_{3}$(t)','$s_{4}$(t)','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
% print -dpdf Time_Vs_TStates_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Time_Vs_TStates_RL_BF_FCL.mat') %saving data in the MAT format

% State trajectories of the estimated value of the parameters
figure(2)
plot(t,z(:,P.n+1:P.n+P.l),'Linewidth',2);
xlabel('Time, t(s)','interpreter','latex')
ylabel('Estimated Theta, s (t)','interpreter','latex')
set(gca,'FontSize',10)
l=legend('$\theta_{1}$(t)','$\theta_{2}$(t)','$\theta_{3}$(t)','$\theta_{4}$(t)','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 


% State trajectories of the estimated value of the Actor weights 
figure(3)
plot(t,z(:,P.n+P.l+P.n*P.l+P.l^2+P.l+1:P.n+P.l+P.n*P.l+P.l^2+P.l+P.W),'Linewidth',2);
xlabel('Time, t(s)','interpreter','latex')
ylabel('Estimated Actor weights, W_{a} (t)','interpreter','latex')
set(gca,'FontSize',10)
l=legend('$W_{a_{1}}$(t)','$W_{a_{2}}$(t)','$W_{a_{3}}$(t)','$W_{a_{4}}$(t)','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 

% State trajectories of the estimated value of the Critic weights 
figure(4)
plot(t,z(:,P.n+P.l+P.n*P.l+P.l^2+P.l+P.W+1:P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W),'Linewidth',2);
xlabel('Time, t(s)','interpreter','latex')
ylabel('Estimated Critic weights, W_{c} (t)','interpreter','latex')
set(gca,'FontSize',10)
l=legend('$W_{c_{1}}$(t)','$W_{c_{2}}$(t)','$W_{c_{3}}$(t)','$W_{c_{4}}$(t)','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 

% State trajectories of the x-coordinates
figure(5)
plot(t,z(:,P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+P.G^2+1:P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+P.G^2+P.n),'Linewidth',2);
xlabel('Time, t(s)','interpreter','latex')
ylabel('Original States, x (t)','interpreter','latex')
set(gca,'FontSize',10)
l=legend('$x_{1}$(t)','$x_{2}$(t)','$x_{3}$(t)','$x_{4}$(t)','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
% print -dpdf Time_Vs_TStates_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Time_Vs_TStates_RL_BF_FCL.mat') %saving data in the MAT format

%checking value of Critic Weight
A = z(end,55:64)
toc

function zDot = closedLoopDynamics(t,z,P)
% Control Penalty Matrix  
R=diag([1,1]); 
% State Penalty Matrix  
Q=diag([1,1,1,1]); 

% Tuning Gains>> Please check the Frontier paper to understand the notation
etac1 = 0.1;
etac2 = 10;
etaa1 = 20; %increase to reduce the difference between Wc and Wa
etaa2 = 0.2;
beta = 0.8;
v1 = 125;


a1 = -7; A1 = 5; %Barrier Boundaries for x1, x1-coordinate values will be remained within (a1,A1) 
a2 = -7; A2 = 5; %Barrier Boundaries for x2, x2-coordinate values will be remained within (a2,A2) 
a3 = -5; A3 = 7; %Barrier Boundaries for x3, x3-coordinate values will be remained within (a3,A3) 
a4 = -5; A4 = 7; %Barrier Boundaries for x4, x4-coordinate values will be remained within (a4,A4) 

% known parameters, 
p1 = 3.473;
p2 = 0.196;
p3 = 0.242;

% Optimizing Matrix/Integraing matrix >> This z matrix is the integration
% of the zDot matrix given at the end of this coding.  
s = z(1:P.n,1);% row 1-4 in z vector
thetaHat = z(P.n+1:P.n+P.l,1);% row 5-8 in z vector 
scriptY = reshape(z(P.n+P.l+1:P.n+P.l+P.n*P.l),P.n,P.l); % row 9-24 in z vector
scriptYf = reshape(z(P.n+P.l+P.n*P.l+1:P.n+P.l+P.n*P.l+P.l^2),P.l,P.l); % row 25-40 in z vector
scriptXf = z(P.n+P.l+P.n*P.l+P.l^2+1:P.n+P.l+P.n*P.l+P.l^2+P.l); % row 41-44 in z vector 
WaH = z(P.n+P.l+P.n*P.l+P.l^2+P.l+1:P.n+P.l+P.n*P.l+P.l^2+P.l+P.W); % row 45-54 in z vector 
WcH = z(P.n+P.l+P.n*P.l+P.l^2+P.l+P.W+1:P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W); % row 55-64 in z vector 
Gamma = reshape(z(P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+1:P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+P.G^2),P.G,P.G); % row 65-164 in z vector 
a = z(P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+P.G^2+1:P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+P.G^2+P.n); % row 165-168 in z vector 
int_G = z(P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+P.G^2+P.n+1:P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+P.G^2+P.n+P.G1); % row 169-172 in z vector 
Cost= z(P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+P.G^2+P.n+P.G1+1); % row 173 in z vector 

%% BE Extrapolation - Terms in weight update laws that are evaluated along arbitrarily selected trajectories xi

%Creating a meshgrid to train to get the optimal actor weight and critic weight along the meshgrid
%points using BE extrapolation

tmp = linspace(-4,3,10);
tmp = meshgrid(tmp,tmp);
tmpt = tmp.';
N = numel(tmp);
S = [tmp(:) tmpt(:) tmpt(:) tmpt(:)];

%initialization of the trained weights in the meshgrid 
clWc=zeros(size(WcH)); %initialization of the critic weights in the meshgrid
clWa=zeros(size(WaH)); %initialization of the actor weights in the meshgrid
CLmatrix=zeros(size(WcH,1),size(WcH,1)); %Matrixized
CLGamma=zeros(size(Gamma)); %initialization of the Gamma weights in the meshgrid


for i=1:N
theta = [5.3;1.1;8.45;2.35]; %unknown parameters
sii = (S(i,:)).';  % Vectorized


a1 = -7; A1 = 5; %Barrier Boundaries for x1, x1-coordinate values will be remained within (a1,A1) 
a2 = -7; A2 = 5; %Barrier Boundaries for x2, x2-coordinate values will be remained within (a2,A2) 
a3 = -5; A3 = 7; %Barrier Boundaries for x3, x3-coordinate values will be remained within (a3,A3) 
a4 = -5; A4 = 7; %Barrier Boundaries for x4, x4-coordinate values will be remained within (a4,A4) 

% known parameters, 
p1 = 3.473;
p2 = 0.196;
p3 = 0.242;

%transforming initial x-coordinates to initial s-coordinates. 
s1 = (log ((A1/a1)*((a1-sii(1))/(A1-sii(1)))));
s2 = (log ((A2/a2)*((a2-sii(2))/(A2-sii(2)))));
s3 = (log ((A3/a3)*((a3-sii(3))/(A3-sii(3)))));
s4 = (log ((A4/a4)*((a4-sii(4))/(A4-sii(4)))));
si  = [s1;s2;s3;s4];

%transforming initial s-coordinates to initial x-coordinates. 
x1_si = a1*A1*((exp(si(1)))-1)/((a1*exp(si(1)))-A1);
x2_si = a2*A2*((exp(si(2)))-1)/((a2*exp(si(2)))-A2);
x3_si = a3*A3*((exp(si(3)))-1)/((a3*exp(si(3)))-A3);
x4_si = a4*A4*((exp(si(4)))-1)/((a4*exp(si(4)))-A4);

%% Defining dynamics in the meshgrid, check ArxiV paper to understand the dynamics
c2 = cos(x2_si);
s2 = sin(x2_si);
M = [p1+2*p3*c2, p2+p3*c2;p2+p3*c2, p2 ];
V = [-p3*s2*x4_si, -p3*s2*(x3_si+x4_si); p3*s2*x3_si, 0 ];
L = inv(M);
v5 = [x3_si,x4_si,tanh(x3_si),tanh(x4_si)];
D = diag(v5);
ll = [          -(p2*x3_si)/(c2^2*p3^2 + p2^2 - p1*p2),    (x4_si*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2), ...
    -(p2*tanh(x3_si))/(c2^2*p3^2 + p2^2 - p1*p2),    (tanh(x4_si)*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2)];
gg = [ (x3_si*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2), -(x4_si*(p1 + 2*c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2), ...
    (tanh(x3_si)*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2), -(tanh(x4_si)*(p1 + 2*c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2)];

hh =  - x3_si*((p3*s2*x3_si*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2)...
    + (p2*p3*s2*x4_si)/(c2^2*p3^2 + p2^2 - p1*p2)) - (p2*p3*s2*x4_si*(x3_si + x4_si))/(c2^2*p3^2 + p2^2 - p1*p2);
jj = x3_si*((p3*s2*x3_si*(p1 + 2*c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2) ...
    + (p3*s2*x4_si*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2)) ...
    + (p3*s2*x4_si*(p2 + c2*p3)*(x3_si + x4_si))/(c2^2*p3^2 + p2^2 - p1*p2);

%Note: We do not need to do the barrier transformation in the meshgrid to train
%our controller as this method is an extrapolation method, and we are picking some random
%points (around the meshgrid) to train. Therefore, we just need to create a
%meshgrid to train, but here we are doing the barrier transformation just
%to mimic the environment, not sure whether it gives a better result or
%not. 


%Dynamics in the meshgrid 
F1_1_si =  (((a1^2*exp(si(1))) - (2*a1*A1) + (A1^2 * exp (-si(1)))) / (-a1*A1^2+A1*a1^2)) * x3_si;
F2_1_si=  (((a2^2*exp(si(2))) - (2*a2*A2) + (A2^2 * exp (-si(2)))) / ( A2*a2^2 -a2*A2^2)) * x4_si;
F3_1_si = (((a3^2*exp(si(3))) - (2*a3*A3) + (A3^2 * exp (-si(3)))) / ( A3*a3^2 -a3*A3^2))* hh;
F4_1_si = (((a4^2*exp(si(4))) - (2*a4*A4) + (A4^2 * exp (-si(4)))) / ( A4*a4^2 -a4*A4^2))* jj;
F3_2_si = -(((a3^2*exp(si(3))) - (2*a3*A3) + (A3^2 * exp (-si(3)))) / ( A3*a3^2 -a3*A3^2))* ll;
F4_2_si = -(((a4^2*exp(si(4))) - (2*a4*A4) + (A4^2 * exp (-si(4)))) / ( A4*a4^2 -a4*A4^2))* gg;

g_3_si = (((a3^2*exp(si(3))) - (2*a3*A3) + (A3^2 * exp (-si(3)))) / ( A3*a3^2 -a3*A3^2))*...
   [ -p2/(c2^2*p3^2 + p2^2 - p1*p2),    (p2 + c2*p3)/(c2^2*p3^2 + p2^2 - p1*p2)]; 
g_4_si = (((a4^2*exp(si(4))) - (2*a4*A4) + (A4^2 * exp (-si(4)))) / ( A4*a4^2 -a4*A4^2))*...
    [ (p2 + c2*p3)/(c2^2*p3^2 + p2^2 - p1*p2), -(p1 + 2*c2*p3)/(c2^2*p3^2 + p2^2 - p1*p2)]; 
Y1_si = [F1_1_si;F2_1_si;F3_1_si;F4_1_si];
Y2_si = [0,0,0,0;0,0,0,0;F3_2_si;F4_2_si];
  G_si = [0,0;0,0;g_3_si;g_4_si];

%ADP-Controller
% Jacobian of the basis vector in the meshgrid
phi_p_si = [si(3),0,si(1),0;0,si(4),0,si(2);...
    0,si(3),si(2),0;si(4),0,0,si(1);...
    si(2),si(1),0,0;0,0,si(4),si(3);...
    2*si(1),0,0,0;0,2*si(2),0,0;...
    0,0,2*si(3),0;0,0,0,2*si(4)];

mu_si = -0.5.*(R\G_si')*phi_p_si'*WaH; % trained controller of the meshgrid

% Cost function in the meshgrid 
ri = si'*Q*si + mu_si'*R*mu_si;

%Definition
omegai = phi_p_si*(Y1_si+Y2_si*thetaHat+G_si*mu_si);
Gsigmai = phi_p_si*G_si*(R\G_si')*phi_p_si';
deltai = WcH'*omegai + ri;
rhoi = (1+v1*(omegai'*omegai));

% actor weight, critic weight, Gamma update law  
clWc = clWc + Gamma*omegai*deltai/rhoi;
clWa = clWa + Gsigmai'*WaH*omegai'*WcH/(4*rhoi);
CLmatrix = CLmatrix + omegai*omegai'/(rhoi);
CLGamma = CLGamma + Gamma*(omegai*omegai'/rhoi^2)*Gamma;
end
clWc=-etac2*clWc/N;
clWa=etac2*clWa/N;
CLGamma=etac2*CLGamma/N;
cbar = min(svd(CLmatrix));


%% Dynamics Section

%transforming initial s-coordinates to initial x-coordinates. 
x1 = a1*A1*((exp(s(1)))-1)/((a1*exp(s(1)))-A1);
x2 = a2*A2*((exp(s(2)))-1)/((a2*exp(s(2)))-A2);
x3 = a3*A3*((exp(s(3)))-1)/((a3*exp(s(3)))-A3);
x4 = a4*A4*((exp(s(4)))-1)/((a4*exp(s(4)))-A4);

%% Defining dynamics, check ArxiV paper to understand the dynamics
c2 = cos(x2);
s2 = sin(x2);
M = [p1+2*p3*c2, p2+p3*c2;p2+p3*c2, p2 ];
V = [-p3*s2*x4, -p3*s2*(x3+x4); p3*s2*x3, 0 ];
L = inv(M);
v5 = [x3,x4,tanh(x3),tanh(x4)];
D = diag(v5);
ll = [          -(p2*x3)/(c2^2*p3^2 + p2^2 - p1*p2),    (x4*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2), ...
    -(p2*tanh(x3))/(c2^2*p3^2 + p2^2 - p1*p2),    (tanh(x4)*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2)];

gg = [ (x3*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2), -(x4*(p1 + 2*c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2), ...
    (tanh(x3)*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2), -(tanh(x4)*(p1 + 2*c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2)];


hh =  - x3*((p3*s2*x3*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2)...
    + (p2*p3*s2*x4)/(c2^2*p3^2 + p2^2 - p1*p2)) - (p2*p3*s2*x4*(x3 + x4))/(c2^2*p3^2 + p2^2 - p1*p2);
jj = x3*((p3*s2*x3*(p1 + 2*c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2) ...
    + (p3*s2*x4*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2)) ...
    + (p3*s2*x4*(p2 + c2*p3)*(x3 + x4))/(c2^2*p3^2 + p2^2 - p1*p2);

%Dynamics
F1_1 = (((a1^2*exp(s(1))) - (2*a1*A1) + (A1^2 * exp (-s(1)))) / (A1*a1^2-a1*A1^2)) * x3;
F2_1 = (((a2^2*exp(s(2))) - (2*a2*A2) + (A2^2 * exp (-s(2)))) / ( A2*a2^2 -a2*A2^2)) * x4;
F3_1 = (((a3^2*exp(s(3))) - (2*a3*A3) + (A3^2 * exp (-s(3)))) / ( A3*a3^2 -a3*A3^2))* hh;
F4_1 = (((a4^2*exp(s(4))) - (2*a4*A4) + (A4^2 * exp (-s(4)))) / ( A4*a4^2 -a4*A4^2))* jj;

F3_2 = -(((a3^2*exp(s(3))) - (2*a3*A3) + (A3^2 * exp (-s(3)))) / ( A3*a3^2 -a3*A3^2))* ll;
F4_2 = -(((a4^2*exp(s(4))) - (2*a4*A4) + (A4^2 * exp (-s(4)))) / ( A4*a4^2 -a4*A4^2))* gg;

g_3 = (((a3^2*exp(s(3))) - (2*a3*A3) + (A3^2 * exp (-s(3)))) / ( A3*a3^2 -a3*A3^2))*...
    [          -p2/(c2^2*p3^2 + p2^2 - p1*p2),    (p2 + c2*p3)/(c2^2*p3^2 + p2^2 - p1*p2)]; 
g_4 = (((a4^2*exp(s(4))) - (2*a4*A4) + (A4^2 * exp (-s(4)))) / ( A4*a4^2 -a4*A4^2))*...
    [ (p2 + c2*p3)/(c2^2*p3^2 + p2^2 - p1*p2), -(p1 + 2*c2*p3)/(c2^2*p3^2 + p2^2 - p1*p2)]; 



%% Defining again to shorten the notations
F1_11  =  x3;
F2_11  =  x4;
F3_11  =  hh;
F4_11  =  jj;
F3_21 =  -ll;
F4_21 = -gg;
g_31 =     [          -p2/(c2^2*p3^2 + p2^2 - p1*p2),    (p2 + c2*p3)/(c2^2*p3^2 + p2^2 - p1*p2)]; 
g_41 =     [ (p2 + c2*p3)/(c2^2*p3^2 + p2^2 - p1*p2), -(p1 + 2*c2*p3)/(c2^2*p3^2 + p2^2 - p1*p2)]; 

%% After Defining all the dynamics we get,  
%% Transformed Dynamics
Y1 = [F1_1;F2_1;F3_1;F4_1];
Y2 = [0,0,0,0;0,0,0,0;F3_2;F4_2];
G1 = [0,0;0,0;g_3;g_4];

%% Original Dyanmics
Y3 = [F1_11;F2_11;F3_11;F4_11];
Y4 = [0,0,0,0;0,0,0,0;F3_21;F4_21];
G2 = [0,0;0,0;g_31;g_41];
theta = [5.3;1.1;8.45;2.35];

%% Terms in weight update laws that are evaluated along the system trajectory x

% Jacobian of basis vector
phi_p = [s(3),0,s(1),0;0,s(4),0,s(2);...
    0,s(3),s(2),0;s(4),0,0,s(1);...
    s(2),s(1),0,0;0,0,s(4),s(3);...
    2*s(1),0,0,0;0,2*s(2),0,0;...
    0,0,2*s(3),0;0,0,0,2*s(4)];

% controller of the transformed system
mu = -0.5.*(R\G1')*phi_p'*WaH;

%Cost Function in the s-coordinate
r = s'*Q*s + mu'*R*mu;

%Definition
omega = phi_p*(Y1+Y2*thetaHat+G1*mu);
Gsigma = phi_p*G1*(R\G1')*phi_p';
delta = WcH'*omega + r;
rho = (1+v1*(omega'*omega));



maxnorm = 100; % Threshold to turn off filters
update = 1;

% Checking norm condition of the scriptYf: 
if norm(scriptYf) > maxnorm
    update = 0;
end

% tuning maxtrix 
GammaM= [100,0,0,0;0,100,0,0;0,0,100,0;0,0,0,100];

%ODE function
zDot = [Y1+Y2*theta+G1*mu;
        GammaM*(scriptXf - scriptYf*thetaHat);
       reshape(Y2, P.n*P.l,1)* update;
       reshape(scriptY.'*scriptY , P.l^2,1)* update;
    scriptY.'*(s-P.s0-int_G) * update;
      (etaa1*(WcH-WaH)-etaa2*WaH+etac1*Gsigma'*WaH*omega'*WcH/(4*rho)+clWa)*update;
      (-etac1*Gamma*omega*delta/rho+clWc)*update;
      reshape(beta*Gamma-CLGamma-etac1*Gamma*(omega*omega'/rho^2)*Gamma,P.G^2,1)*update;
      Y3+Y4*theta+G2*mu;
      Y1+G1*mu;
      r;];
end

