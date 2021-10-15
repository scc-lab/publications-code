%% Paper title "Safe Model-Based Reinforcement Learning for Systems with Parametric Uncertainties" 
%% arXiv:2007.12666
%% Coded by S M Nahid Mahmud, MS Grad Student, Oklahoma State University.
%% nahid.mahmud@okstate.edu
%% Four-state dynamical system (Section 7.2 of the journal paper)

function phaseout = continuous(input)
a1 = -7; A1 = 5; %Barrier Boundaries for x1, x1-coordinate values will be remained within (a1,A1) 
a2 = -7; A2 = 5; %Barrier Boundaries for x2, x2-coordinate values will be remained within (a2,A2) 
a3 = -5; A3 = 7; %Barrier Boundaries for x3, x3-coordinate values will be remained within (a3,A3) 
a4 = -5; A4 = 7; %Barrier Boundaries for x4, x4-coordinate values will be remained within (a4,A4) 

% known parameters, 
p1 = 3.473;
p2 = 0.196;
p3 = 0.242;

t = input.phase.time; % Denoting symbol for input time 
x = input.phase.state; % Denoting symbol for input state  
u = input.phase.control; % Denoting symbol for input control 

%transforming initial s-coordinates to initial x-coordinates. 
x1_now = a1*A1*((exp(x(:,1)))-1)./((a1*exp(x(:,1)))-A1);
x2_now = a2*A2*((exp(x(:,2)))-1)./((a2*exp(x(:,2)))-A2);
x3_now = a3*A3*((exp(x(:,3)))-1)./((a3*exp(x(:,3)))-A3);
x4_now = a4*A4*((exp(x(:,4)))-1)./((a4*exp(x(:,4)))-A4);


%% Defining dynamics, check ACC paper to understand the dynamics
c2 = cos(x2_now);
s2 = sin(x2_now);

hh =  - x3_now.*((p3.*s2.*x3_now.*(p2 + c2*p3))./(c2.^2*p3^2 + p2^2 - p1*p2)...
    + (p2*p3.*s2.*x4_now)./(c2.^2*p3^2 + p2^2 - p1*p2)) - (p2*p3*s2.*x4_now.*(x3_now + x4_now))./(c2.^2*p3^2 + p2^2 - p1*p2);

jj = x3_now.*((p3.*s2.*x3_now.*(p1 + 2.*c2.*p3))./(c2.^2*p3^2 + p2^2 - p1*p2) ...
    + (p3*s2.*x4_now.*(p2 + c2.*p3))./(c2.^2.*p3^2 + p2^2 - p1*p2)) ...
    + (p3.*s2.*x4_now.*(p2 + c2.*p3).*(x3_now + x4_now))./(c2.^2.*p3^2 + p2^2 - p1*p2);

ll_1 =  -((p2.*x3_now)./(c2.^2.*p3^2 + p2^2 - p1*p2))*5.3 + ((x4_now.*(p2 + c2.*p3))./(c2.^2.*p3^2 + p2^2 - p1*p2))*1.1  -((p2.*tanh(x3_now))./(c2.^2.*p3^2 + p2^2 - p1*p2))*8.45 ...
    + ((tanh(x4_now).*(p2 + c2.*p3))./(c2.^2.*p3^2 + p2^2 - p1*p2))*2.35;

gg_1 = ((x3_now.*(p2 + c2.*p3))./(c2.^2.*p3^2 + p2^2 - p1*p2))*5.3 - ((x4_now.*(p1 + 2.*c2.*p3))./(c2.^2.*p3^2 + p2^2 - p1*p2))*1.1 + ((tanh(x3_now).*(p2 + c2.*p3)) ./(c2.^2.*p3^2 + p2^2 - p1*p2))*8.45...
         - ((tanh(x4_now).*(p1 + 2.*c2.*p3))./(c2.^2.*p3^2 + p2^2 - p1*p2))*2.35;
     
ee_1  = (-p2./(c2.^2.*p3^2 + p2^2 - p1*p2)).*(u(:,1)) + ((p2 + c2.*p3)./(c2.^2.*p3^2 + p2^2 - p1*p2)).*(u(:,2)) ; 

ff_1 = ((p2 + c2.*p3)./(c2.^2.*p3^2 + p2^2 - p1*p2)).*(u(:,1)) + (-(p1 + 2.*c2.*p3)./(c2.^2.*p3^2 + p2^2 - p1*p2)).*(u(:,2));


%% Transformed Dynamics
F1_1 = ((A1.^2*exp(-x(:,1))-2*a1*A1+a1.^2*exp(x(:,1))) ./ ( A1*a1.^2 -a1*A1.^2)) .* x3_now;
F2_1 = ((A2.^2*exp(-x(:,2))-2*a2*A1+a2.^2*exp(x(:,2))) ./ ( A2*a2.^2 -a2*A2.^2)) .* x4_now;
F3_1 = ((A3.^2*exp(-x(:,3))-2*a3*A1+a3.^2*exp(x(:,3))) ./ ( A3*a3.^2 -a3*A3.^2)) .* hh;
F4_1 = ((A4.^2*exp(-x(:,4))-2*a4*A1+a4.^2*exp(x(:,4))) ./ ( A4*a4.^2 -a4*A4.^2)) .* jj;
F3_2 = -((A3.^2*exp(-x(:,3))-2*a3*A1+a3.^2*exp(x(:,3))) ./ ( A3*a3.^2 -a3*A3.^2)) .*  ll_1;
F4_2 = -((A4.^2*exp(-x(:,4))-2*a4*A1+a4.^2*exp(x(:,4))) ./ ( A4*a4.^2 -a4*A4.^2)) .* gg_1;
g_3 = ((A3.^2*exp(-x(:,3))-2*a3*A1+a3.^2*exp(x(:,3))) ./ ( A3*a3.^2 -a3*A3.^2)) .* ee_1;
g_4 = ((A4.^2*exp(-x(:,4))-2*a4*A1+a4.^2*exp(x(:,4))) ./ ( A4*a4.^2 -a4*A4.^2)) .* ff_1;
    

%% ODE function
xdot = [F1_1, F2_1, F3_1+F3_2+g_3, F4_1+F4_2+g_4];
phaseout.dynamics = xdot;
phaseout.integrand = x(:,1).^2 + x(:,2).^2 + x(:,3).^2 + x(:,4).^2 + (u(:,1).^2)+(u(:,2).^2); %Cost function