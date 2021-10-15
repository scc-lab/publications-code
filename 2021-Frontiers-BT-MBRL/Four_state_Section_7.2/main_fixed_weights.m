%% Paper title "Safe Model-Based Reinforcement Learning for Systems with Parametric Uncertainties" 
%% arXiv:2007.12666
%% Coded by S M Nahid Mahmud, MS Grad Student, Oklahoma State University.
%% nahid.mahmud@okstate.edu
%% Four-state dynamical system (Section 7.2 of the journal paper)

clear all
clc
tic

%% PreSim
P.n = 4; %number of states
P.m = 1; %number of dimension
P.l = 4; %number of unknown parameters
P.x0 = [-5;-5;5;5]; % Initial x-coordinate values, it has to be within (a,A).

% known parameters, 
p1 = 3.473;
p2 = 0.196;
p3 = 0.242;


a1 = -7; A1 = 5; %Barrier Boundaries for x1, x1-coordinate values will be remained within (a1,A1) 
a2 = -7; A2 = 5; %Barrier Boundaries for x2, x2-coordinate values will be remained within (a2,A2) 
a3 = -5; A3 = 7; %Barrier Boundaries for x3, x3-coordinate values will be remained within (a3,A3) 
a4 = -5; A4 = 7; %Barrier Boundaries for x4, x4-coordinate values will be remained within (a4,A4) 


%transforming initial x-coordinates to initial s-coordinates. 
s10 = (log ((A1/a1)*((a1-P.x0(1))/(A1-P.x0(1)))));
s20 = (log ((A2/a2)*((a2-P.x0(2))/(A2-P.x0(2)))));
s30 = (log ((A3/a3)*((a3-P.x0(3))/(A3-P.x0(3)))));
s40 = (log ((A4/a4)*((a4-P.x0(4))/(A4-P.x0(4)))));

P.s0 = [s10;s20;s30;s40]; % Initial s-coordinate values

P.Cost0=0; %Initial Cost 

%Initial values for x = 1-4, Cost = 5
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
l=legend('$s_{1}$(t)','$s_{2}$(t)','$s_{3}$(t)','$s_{4}$(t)','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 

% print -dpng States_RL_BF_robot22.pdf
% save('States_RL_BF_manipulator.mat','-v7.3') %saving data 



function zDot = closedLoopDynamics(t,z,P)
% Control Penalty Matrix  
R=diag([1,1]); 
% State Penalty Matrix  
Q=diag([1,1,1,1]); 


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
s = z(1:P.n,1); % row 1-4 in z vector
Cost= z(P.n+1) % row 5 in z vector

%% Dynamics

%transforming initial s-coordinates to initial x-coordinates. 
x1 = a1*A1*((exp(s(1)))-1)/((a1*exp(s(1)))-A1);
x2 = a2*A2*((exp(s(2)))-1)/((a2*exp(s(2)))-A2);
x3 = a3*A3*((exp(s(3)))-1)/((a3*exp(s(3)))-A3);
x4 = a4*A4*((exp(s(4)))-1)/((a4*exp(s(4)))-A4);
%% Defining dynamics, check ACC paper to understand the dynamics

c2 = cos(x2);
s2 = sin(x2);
M = [p1+2*p3*c2, p2+p3*c2;p2+p3*c2, p2 ];
V = [-p3*s2*x4, -p3*s2*(x3+x4); p3*s2*x3, 0 ];
L = inv(M);
v5 = [x3,x4,tanh(x3),tanh(x4)];
D = diag(v5);
% 
ll = [          -(p2*x3)/(c2^2*p3^2 + p2^2 - p1*p2),    (x4*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2), ...
    -(p2*tanh(x3))/(c2^2*p3^2 + p2^2 - p1*p2),    (tanh(x4)*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2)];

gg = [ (x3*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2), -(x4*(p1 + 2*c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2), ...
    (tanh(x3)*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2), -(tanh(x4)*(p1 + 2*c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2)];


hh =  - x3*((p3*s2*x3*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2)...
    + (p2*p3*s2*x4)/(c2^2*p3^2 + p2^2 - p1*p2)) - (p2*p3*s2*x4*(x3 + x4))/(c2^2*p3^2 + p2^2 - p1*p2);
jj = x3*((p3*s2*x3*(p1 + 2*c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2) ...
    + (p3*s2*x4*(p2 + c2*p3))/(c2^2*p3^2 + p2^2 - p1*p2)) ...
    + (p3*s2*x4*(p2 + c2*p3)*(x3 + x4))/(c2^2*p3^2 + p2^2 - p1*p2);



%% Defining again to shorten the notations

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



%% Transformed Dynamics
Y1 = [F1_1;F2_1;F3_1;F4_1];
Y2 = [0,0,0,0;0,0,0,0;F3_2;F4_2];
G1 = [0,0;0,0;g_3;g_4];

theta = [5.3;1.1;8.45;2.35];



%% Terms in weight update laws that are evaluated along the system trajectory x

% Jacobian of basis vector
phi_p = [s(3),0,s(1),0;0,s(4),0,s(2);...
    0,s(3),s(2),0;s(4),0,0,s(1);...
    s(2),s(1),0,0;0,0,s(4),s(3);...
    2*s(1),0,0,0;0,2*s(2),0,0;...
    0,0,2*s(3),0;0,0,0,2*s(4)];

% Fixed Weights
W = [50.7057    2.1208    0.5059    2.6475    5.7407   -0.3492   45.3566    4.8677   -1.6964    1.1666]';

% controller of the transformed system
mu = -0.5.*(R\G1')*phi_p'*W;

%Cost Function in the s-coordinate
r = s'*Q*s + mu'*R*mu;


%ODE function
zDot = [Y1+Y2*theta+G1*mu;
    r];
end

