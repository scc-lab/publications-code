%% Paper title "Safe Model-Based Reinforcement Learning for Systems with Parametric Uncertainties" 
%% arXiv:2007.12666
%% Coded by S M Nahid Mahmud, MS Grad Student, Oklahoma State University.
%% nahid.mahmud@okstate.edu
%% Two-state dynamical system (Section 7.1 of the journal paper)

function phaseout = continuous(input)
a1 = -7; A1 = 5; %Barrier Boundaries for x1, x1-coordinate values will be remained within (a1,A1) 
a2 = -5; A2 = 7; %Barrier Boundaries for x2, x2-coordinate values will be remained within (a2,A2) 


t = input.phase.time; % Denoting symbol for input time 
x = input.phase.state; % Denoting symbol for input state  
u = input.phase.control; % Denoting symbol for input control 


%transforming initial s-coordinates to initial x-coordinates.
x1_now = a1*A1*((exp(x(:,1)))-1)./((a1*exp(x(:,1)))-A1);
x2_now = a2*A2*((exp(x(:,2)))-1)./((a2*exp(x(:,2)))-A2);


%% Transformed Dynamics, check ACC paper to understand the dynamics

F1_2 =  ((A1.^2*exp(-x(:,1))-2*a1*A1+a1.^2*exp(x(:,1))) ./ ( A1*a1.^2 -a1*A1.^2)) .* x2_now;
f =  - x1_now - 0.5 .* x2_now + 0.5 .* (x2_now.*(cos(2.*x1_now)+2).^2); 
g = (cos(2.*x1_now)+2);
F2 = ((A2.^2 *exp(-x(:,2))-2*a2*A2+a2.^2 *exp(x(:,2))) ./ ( A2*a2.^2 -a2*A2.^2)) .* (f + g .* u );

%% ODE function
xdot = [F1_2,F2];
phaseout.dynamics = xdot;
phaseout.integrand = x(:,1).^2 + x(:,2).^2+u.^2;