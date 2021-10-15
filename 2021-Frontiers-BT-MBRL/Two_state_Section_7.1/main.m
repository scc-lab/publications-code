%% Paper title "Safe Model-Based Reinforcement Learning for Systems with Parametric Uncertainties" 
%% arXiv:2007.12666
%% Coded by S M Nahid Mahmud, MS Grad Student, Oklahoma State University.
%% nahid.mahmud@okstate.edu
%% Two-state dynamical system (Section 7.1 of the journal paper)

clear all
clc
tic
%% PreSim

a1 = -7; A1 = 5; %Barrier Boundaries for x1, x1-coordinate values will be remained within (a1,A1) 
a2 = -5; A2 = 7; %Barrier Boundaries for x2: x2-coordinate values will be remained within (a2,A2) 

P.n = 2; %number of states
P.m = 1; %number of dimension
P.l = 4; %number of unknown parameters
P.x0 = [-6.5;6.5]; % Initial x-coordinate values, it has to be within (a,A).

%transforming initial x-coordinates to initial s-coordinates. 
s10 = (log ((A1/a1)*((a1-P.x0(1))/(A1-P.x0(1)))));
s20 = (log ((A2/a2)*((a2-P.x0(2))/(A2-P.x0(2)))));

P.s0 = [s10;s20]; % Initial s-coordinate values
P.W = 3; % number of the Actor Weight vector elements as wellas Critic Weight vector elements.
P.G =3;  % number of the Gamma vector elements. 

P.theta = [1;-1;-0.5;0.5]; % True Value of the unknown parameters 
P.theta0= [0;0;0;0]; % Picked initial value of the unknown parameters set by the user. 
P.theta_tilde0 = P.theta-P.theta0; % Errors between the value of the unknown parameters and initial picked up value of the unknown parameters. 
P.IG = [0;0]; % Initial value of the integrated G function.  


%%Initial values for s = 1-2, thetahat= 3-6, scriptY = 7-14, scriptYf =
%15-30, scriptXf = 31-34, WaH = 35:37, WcH = 38:40, Gamma =  41:49, 
% x = 50:51, Int_G = 52:53
% ThetaHat denotes Estimated value of the unknown parameters, WaH denotes Estimated Actor weight, WcH denotes Estimated Critic weight,
% Int_G denotes the integration of G function. 
P.z0 = [P.s0;P.theta0;zeros(P.n*P.l+P.l^2+P.l,1);0.5*ones(2*P.W,1);1;0;0;0;1;0;0;0;1;P.x0;P.IG];

%% Integration
%ODE45 function
options = odeset('OutputFcn',@odeplot,'OutputSel',1:P.n);
[t,z] = ode45(@(t,z) closedLoopDynamics(t,z,P),[0 10],P.z0,options);



%% Plot

% Phase portrait of the s-coordinates
figure(1)
plot(z(:,1),z(:,2),'-s','Linewidth',2)
xlabel('$s_{1}$(t)','interpreter','latex')
ylabel('$s_{2}$(t)','interpreter','latex')
set(gca,'FontSize',10)
l=legend('Transformed State Trajectory','interpreter','latex');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
% print -dpdf  States_RL_BF_FCL.pdf %saving plot as a pdf 
% save('States_RL_BF_FCL.mat') %saving data in the MAT format

% State trajectories of the s-coordinates
figure(2)
p=plot(t,z(:,1),t,z(:,2),'Linewidth',2);
mrk1={'s','v','o','*','^','x'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time, t(s)','interpreter','latex')
ylabel('Transformed States, s (t)','interpreter','latex')
set(gca,'FontSize',10)
l=legend('$s_{1}$(t)','$s_{2}$(t)','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
% print -dpdf Time_Vs_TStates_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Time_Vs_TStates_RL_BF_FCL.mat') %saving data in the MAT format

% State trajectories of the x1-coordinates
figure(3)
p=plot(t,z(:,50),'-s','LineWidth',1);
xlabel('Time (s)','FontSize',16)
ylabel('$x_{1}$(t)','interpreter','latex')
set(gca,'FontSize',16)
hold on
p1=plot(t,a1*ones(size(t)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(t,A1*ones(size(t)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
l=legend('$x_{1}$ Trajectory','interpreter','latex');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
% print -dpdf  Time_Vs_Barrier_x1_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Time_Vs_Barrier_x1_RL_BF_FCL.mat') %saving data in the MAT format

% State trajectories of the x2-coordinates
figure(4)
p=plot(t,z(:,51),'-s','LineWidth',1);
xlabel('Time (s)','FontSize',16)
ylabel('$x_{2}$(t)','interpreter','latex')
set(gca,'FontSize',16)
hold on
p1=plot(t,a1*ones(size(t)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(t,A2*ones(size(t)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
l=legend('$x_{1}$ Trajectory','interpreter','latex');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
% print -dpdf  Time_Vs_Barrier_x2_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Time_Vs_Barrier_x2_RL_BF_FCL.mat') %saving data in the MAT format


% Trajectory of the Actor and Critic weights 
figure(5)
p=plot(t,z(:,35),t,z(:,36),t,z(:,37),t,z(:,38),t,z(:,39),t,z(:,40),'Linewidth',2);
mrk1={'s','v','o','*','^','x'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time, t(s)','interpreter','latex')
ylabel('Weight Estimations for States, $\hat{W}$','interpreter','latex')
set(gca,'FontSize',10)
l=legend('$\hat{W_{a}}_{1}$','$\hat{W_{a}}_{2}$','$\hat{W_{a}}_{3}$','$\hat{W_{c}}_{1}$','$\hat{W_{c}}_{2}$','$\hat{W_{c}}_{3}$','interpreter','latex');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
% print -dpdf  Actor_Weights_Estimations_for_States_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Actor_Weights_Estimations_for_States_RL_BF_FCL.mat') %saving data in the MAT format
 
% Phase portrait of the x-coordinates
figure(6)
plot(z(:,50),z(:,51),'-s','Linewidth',2)
rectangle('position',[-7 -5 12 12])
axis([-8 7 -7 8])
xlabel('$x_{1}$(t)','interpreter','latex')
ylabel('$x_{2}$(t)','interpreter','latex')
set(gca,'FontSize',10)
l=legend('Original State Trajectory','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
% print -dpdf  Main_States_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Main_States_RL_BF_FCL.mat') %saving data in the MAT format



% State trajectories of the x-coordinates
figure(7)
p=plot(t,z(:,50),t,z(:,51),'-s','Linewidth',2);
mrk1={'s','v','o','*','^','x'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time, t(s)','interpreter','latex')
ylabel('Original States, s(t)','interpreter','latex')
set(gca,'FontSize',10)
l=legend('$x_{1}$(t)','$x_{2}$(t)','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
% print -dpdf  Time_Vs_OStates_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Time_Vs_OStates_RL_BF_FCL.mat') %saving data in the MAT format

% State trajectories of the parameter estimation errors
figure(8)
p=plot(t,(P.theta(1)-z(:,3)),t,(P.theta(2)-z(:,4)),t,(P.theta(3)-z(:,5)),t,(P.theta(4)-z(:,6)),'Linewidth',2);
mrk1={'s','v','o','*','^','x'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time, t(s)','interpreter','latex')
ylabel('Parameter Estimation Errors, $\tilde{\theta}$','interpreter','latex')
set(gca,'FontSize',10)
l=legend('$\tilde{\theta}_{1}$','$\tilde{\theta}_{2}$','$\tilde{\theta}_{3}$','$\tilde{\theta}_{4}$','interpreter','latex');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
% print -dpdf  Parameter_Estimation_Errors_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Parameter_Estimation_Errors_RL_BF_FCL.mat') %saving data in the MAT format



function zDot = closedLoopDynamics(t,z,P)
% Control Penalty Matrix  
R=0.1; 
% State Penalty Matrix  
Q=[10 0; 0 10]; 
% Tuning Gains
etac1 = 0.3;
etac2 = 5;
etaa1 = 180;
etaa2 = 0.0001;
beta = 0.03;
v = 0.5;

% Optimizing Matrix/Integraing matrix >> This z matrix is the integration
% of the zDot matrix given at the end of this coding.  
s = z(1:P.n,1); % row 1-2 in z vector
thetaHat = z(P.n+1:P.n+P.l,1); % row 3-6 in z vector
scriptY = reshape(z(P.n+P.l+1:P.n+P.l+P.n*P.l),P.n,P.l); % row 7-14 in z vector
scriptYf = reshape(z(P.n+P.l+P.n*P.l+1:P.n+P.l+P.n*P.l+P.l^2),P.l,P.l); % row 15-30 in z vector
scriptXf = z(P.n+P.l+P.n*P.l+P.l^2+1:P.n+P.l+P.n*P.l+P.l^2+P.l); % row 31-34 in z vector
WaH = z(P.n+P.l+P.n*P.l+P.l^2+P.l+1:P.n+P.l+P.n*P.l+P.l^2+P.l+P.W); % row 35-37 in z vector
WcH = z(P.n+P.l+P.n*P.l+P.l^2+P.l+P.W+1:P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W); % row 38-40 in z vector
Gamma = reshape(z(P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+1:P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+P.G^2),P.G,P.G); % row 41-49 in z vector
a = z(P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+P.G^2+1:P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+P.G^2+P.n); % row 50-51 in z vector
int_G = z(P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+P.G^2+P.n+1:P.n+P.l+P.n*P.l+P.l^2+P.l+2*P.W+P.G^2+P.n+2); % row 52-53 in z vector
 
%% BE Extrapolation - Terms in weight update laws that are evaluated along arbitrarily selected trajectories sii

%Creating a meshgrid to train to get the optimal actor weight and critic weight along the meshgrid
%points using BE extrapolation
tmp = linspace(-4,3,10);
tmp = meshgrid(tmp,tmp);
tmpt = tmp.';
M = numel(tmp);
S = [tmp(:) tmpt(:)];

%initialization of the trained weights in the meshgrid 
clWc=zeros(size(WcH)); %initialization of the critic weights in the meshgrid
clWa=zeros(size(WaH)); %initialization of the actor weights in the meshgrid
CLmatrix=zeros(size(WcH,1),size(WcH,1)); %Matrixized
CLGamma=zeros(size(Gamma)); %initialization of the Gamma weights in the meshgrid


for i=1:M
sii = (S(i,:)).'; % Vectorized

a1 = -7; A1 = 5; %Barrier Boundaries for sii(1), sii(1)-coordinate values will be remained within (a1,A1) 
a2 = -5; A2 = 7; %Barrier Boundaries for sii(2), sii(2)-coordinate values will be remained within (a2,A2) 

%transforming initial x-coordinates to initial s-coordinates in the meshgrid. 
s1 = (log ((A1/a1)*((a1-sii(1))/(A1-sii(1)))));
s2 = (log ((A2/a2)*((a2-sii(2))/(A2-sii(2)))));
si  = [s1;s2];

%transforming initial s-coordinates to initial x-coordinates in the meshgrid. 
x1_si = a1*A1*((exp(si(1)))-1)/((a1*exp(si(1)))-A1);
x2_si = a2*A2*((exp(si(2)))-1)/((a2*exp(si(2)))-A2);


%Note: We do not need to do the barrier transformation in the meshgrid to train
%our controller as this method is an extrapolation method, and we are picking some random
%points (around the meshgrid) to train. Therefore, we just need to create a
%meshgrid to train, but here we are doing the barrier transformation just
%to mimic the environment, not sure whether it gives a better result or
%not. 


%Dynamics in the meshgrid 
F1_1_si =  (((a1^2*exp(si(1))) - (2*a1*A1) + (A1^2 * exp (-si(1)))) / (-a1*A1^2+A1*a1^2)) * x1_si;
F1_2_si =  (((a1^2*exp(si(1))) - (2*a1*A1) + (A1^2 * exp (-si(1)))) / (-a1*A1^2+A1*a1^2)) * x2_si;
F2_1_si=  (((a2^2*exp(si(2))) - (2*a2*A2) + (A2^2 * exp (-si(2)))) / ( A2*a2^2 -a2*A2^2)) * x1_si;
F2_2_si = (((a2^2*exp(si(2))) - (2*a2*A2) + (A2^2 * exp (-si(2)))) / ( A2*a2^2 -a2*A2^2))* x2_si;
F2_3_si = (((a2^2*exp(si(2))) - (2*a2*A2) + (A2^2 * exp (-si(2)))) / ( A2*a2^2 -a2*A2^2)) * (x2_si*(cos(2*x1_si)+2)^2);
G2_si = (((a2^2*exp(si(2))) - (2*a2*A2) + (A2^2 * exp (-si(2)))) / ( A2*a2^2 -a2*A2^2))* (cos(2*x1_si)+2);
Y_si = [F1_2_si,0,0,0;0,F2_1_si,F2_2_si,F2_3_si];
G_si = [0; G2_si];

%ADP-Controller
% Basis vector in the meshgrid: [si(1)^2, si(1)si(2),si(2)^2]
phi_p_si = [2*si(1) 0; si(2) si(1); 0 2*si(2)]; % Jacobian of the basis vector in the meshgrid
mu_si = -0.5.*(R\G_si')*phi_p_si'*WaH; % trained controller of the meshgrid

% Cost function in the meshgrid 
ri = si'*Q*si + mu_si'*R*mu_si; 

%Definition
omegai = phi_p_si*(Y_si*thetaHat+G_si*mu_si);
Gsigmai = phi_p_si*G_si*(R\G_si')*phi_p_si';
deltai = WcH'*omegai + ri;
rhoi = (1+v*(omegai'*omegai));


% actor weight, critic weight, Gamma update law  
clWc = clWc + Gamma*omegai*deltai/rhoi;
clWa = clWa + Gsigmai'*WaH*omegai'*WcH/(4*rhoi);
CLmatrix = CLmatrix + omegai*omegai'/(rhoi);
CLGamma = CLGamma + Gamma*(omegai*omegai'/rhoi^2)*Gamma;
end
clWc=-etac2*clWc/M;
clWa=etac2*clWa/M;
CLGamma=etac2*CLGamma/M;
cbar = min(svd(CLmatrix));

%% Dynamics

%%transforming initial s-coordinates to initial x-coordinates. 
x1 = a1*A1*((exp(s(1)))-1)/((a1*exp(s(1)))-A1);
x2 = a2*A2*((exp(s(2)))-1)/((a2*exp(s(2)))-A2);

%% Original dynamics without the parameters
Y1 = [x2,0,0,0;0,x1,x2, x2*((cos(2*x1)+2)^2)];
G1 = [0; (cos(2*x1)+2)];
%Known Parameters of the dynamics
theta = [1;-1;-.5;.5];


%%Transformed dynamics
F1_2 =  (((a1^2*exp(s(1))) - (2*a1*A1) + (A1^2 * exp (-s(1)))) / ( A1*a1^2 -a1*A1^2)) * x2;
F2_1 =  (((a2^2*exp(s(2))) - (2*a2*A2) + (A2^2 * exp (-s(2)))) / ( A2*a2^2 -a2*A2^2)) * x1;
F2_2 = (((a2^2*exp(s(2))) - (2*a2*A2) + (A2^2 * exp (-s(2)))) / ( A2*a2^2 -a2*A2^2)) * x2;
F2_3 = (((a2^2*exp(s(2))) - (2*a2*A2) + (A2^2 * exp (-s(2)))) / ( A2*a2^2 -a2*A2^2)) * (x2*(cos(2*x1)+2)^2);
G2 = (((a2^2*exp(s(2))) - (2*a2*A2) + (A2^2 * exp (-s(2)))) / ( A2*a2^2 -a2*A2^2)) * (cos(2*x1)+2);
Y = [F1_2,0,0,0;0,F2_1,F2_2,F2_3];
G = [0; G2];


%% Terms in weight update laws that are evaluated along the system trajectory x

% Basis vector: [s1^2, s1s2,s2^2]
phi_p = [2*s(1) 0; s(2) s(1); 0 2*s(2)]; % Jacobian of basis vector

% controller of the transformed system
mu = -0.5.*(R\G')*phi_p'*WaH;

%Cost Function in the s-coordinate
r = s'*Q*s + mu'*R*mu;

%Definition
omega = phi_p*(Y*thetaHat+G*mu);
Gsigma = phi_p*G*(R\G')*phi_p';
delta = WcH'*omega + r;
rho = (1+v*(omega'*omega));



maxNorm = 1000; % Threshold to turn off filters
update = 1;

% Checking norm condition of the scriptYf: 
if norm(scriptYf) > maxNorm
    update = 0;
end

% tuning maxtrix 
% GammaM= [150,0,0,0;0,150,0,0;0,0,150,0;0,0,0,150]; % it can also be used.
GammaM= 150;

%Integrating function
zDot = [Y*theta+G*mu;
         GammaM*(scriptXf - scriptYf*thetaHat);
       reshape(Y , P.n*P.l,1)* update;
       reshape(scriptY.'*scriptY , P.l^2,1)* update;
    scriptY.'*(s-P.s0-int_G) * update
      (etaa1*(WcH-WaH)-etaa2*WaH+etac1*Gsigma'*WaH*omega'*WcH/(4*rho)+clWa)*update;
      (-etac1*Gamma*omega*delta/rho+clWc)*update;
      reshape(beta*Gamma-CLGamma-etac1*Gamma*(omega*omega'/rho^2)*Gamma,P.G^2,1)*update;
      Y1*theta+G1*mu;
      G*mu];
end
