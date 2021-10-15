%% Paper title "Safe Model-Based Reinforcement Learning for Systems with Parametric Uncertainties" 
%% arXiv:2007.12666
%% Coded by S M Nahid Mahmud, MS Grad Student, Oklahoma State University.
%% nahid.mahmud@okstate.edu
%% Two-state dynamical system (Section 7.1 of the journal paper)

clear all
clc
tic

%% PreSim
P.n = 2; %number of states
P.m = 1; %number of dimension
P.l = 4; %number of unknown parameters
P.x0 = [-6.5;6.5]; % Initial x-coordinate values, it has to be within (a,A).



a1 = -7; A1 = 5; %Barrier Boundaries for x1, x1-coordinate values will be remained within (a1,A1) 
a2 = -5; A2 = 7; %Barrier Boundaries for x2: x2-coordinate values will be remained within (a2,A2) 


%transforming initial x-coordinates to initial s-coordinates. 
s10 = (log ((A1/a1)*((a1-P.x0(1))/(A1-P.x0(1)))));
s20 = (log ((A2/a2)*((a2-P.x0(2))/(A2-P.x0(2)))));

P.s0 = [s10;s20]; % Initial s-coordinate values

P.Cost0=0; %Initial Cost 

%Initial values for x = 1-2, Cost = 3
P.z0 = [P.s0;P.Cost0];

%% Integration
%ODE45 function
% options = odeset('OutputFcn',@odeplot,'OutputSel',P.n+1:P.n+P.l);
[t,z] = ode45(@(t,z) closedLoopDynamics(t,z,P),[0 100],P.z0) ;%ODE45 function


% State trajectories of the s-coordinates
figure(1)
plot(t,z(:,1:P.n),'Linewidth',2);
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

% print -dpng States_RL_BF_robot22.pdf
% save('States_RL_BF.mat','-v7.3') %saving data 



function zDot = closedLoopDynamics(t,z,P)
% Control Penalty Matrix  
R=0.1; 
% State Penalty Matrix  
Q=[10 0; 0 10]; 


a1 = -7; A1 = 5; %Barrier Boundaries for sii(1), sii(1)-coordinate values will be remained within (a1,A1) 
a2 = -5; A2 = 7; %Barrier Boundaries for sii(2), sii(2)-coordinate values will be remained within (a2,A2) 


% Optimizing Matrix/Integraing matrix >> This z matrix is the integration
% of the zDot matrix given at the end of this coding.  
s = z(1:P.n,1); % row 1-2 in z vector
Cost= z(P.n+1) % row 3 in z vector

%% Dynamics

%%transforming initial s-coordinates to initial x-coordinates. 
x1 = a1*A1*((exp(s(1)))-1)/((a1*exp(s(1)))-A1);
x2 = a2*A2*((exp(s(2)))-1)/((a2*exp(s(2)))-A2);
%% Defining dynamics, check ACC paper to understand the dynamics

%%Transformed dynamics
F1_2 =  (((a1^2*exp(s(1))) - (2*a1*A1) + (A1^2 * exp (-s(1)))) / ( A1*a1^2 -a1*A1^2)) * x2;
F2_1 =  (((a2^2*exp(s(2))) - (2*a2*A2) + (A2^2 * exp (-s(2)))) / ( A2*a2^2 -a2*A2^2)) * x1;
F2_2 = (((a2^2*exp(s(2))) - (2*a2*A2) + (A2^2 * exp (-s(2)))) / ( A2*a2^2 -a2*A2^2)) * x2;
F2_3 = (((a2^2*exp(s(2))) - (2*a2*A2) + (A2^2 * exp (-s(2)))) / ( A2*a2^2 -a2*A2^2)) * (x2*(cos(2*x1)+2)^2);
G2 = (((a2^2*exp(s(2))) - (2*a2*A2) + (A2^2 * exp (-s(2)))) / ( A2*a2^2 -a2*A2^2)) * (cos(2*x1)+2);

%% Transformed Dynamics
Y1 = [F1_2,0,0,0;0,F2_1,F2_2,F2_3];
G1 = [0; G2];

theta = [1;-1;-.5;.5];



%% Terms in weight update laws that are evaluated along the system trajectory x

% Basis vector: [s1^2, s1s2,s2^2]
phi_p = [2*s(1) 0; s(2) s(1); 0 2*s(2)]; % Jacobian of basis vector

% Fixed Weights
W = [10.7057    0.1208    1.5059]'; 
%Change according to the value of the Actor-Critic weights learnt from the simulation run.

% controller of the transformed system
mu = -0.5.*(R\G1')*phi_p'*W;

%Cost Function in the s-coordinate
r = s'*Q*s + mu'*R*mu;


%ODE function
zDot = [Y1*theta+G1*mu;
    r];
end

