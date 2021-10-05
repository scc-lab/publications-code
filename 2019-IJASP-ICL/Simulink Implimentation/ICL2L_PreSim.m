clear all
% Initial conditions
P.x0 = [1;-1;1;-1];
P.n = 2;
% Control gains
P.k1 = 1;
P.alpha = 1;

% System ID gains
P.Gamma = 0.004*diag([1 1 1 15 20 67 50]);
P.N = 50; % Number of points in storage
P.k2 = 100/P.N;
P.epsilon = 0;
P.Delta_t = 1; % Integration window

P.ts = 0.001; % Sample time

P.thetaStar = [3.473;.196;.242;5.3;1.1;8.45;2.35]; % True parameters
P.theta0 = ones(size(P.thetaStar)); % Initial guess of parameters

P.L = numel(P.thetaStar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Comment out the rest of this file after the first run. %%%%%%%%%%

% Create matlab files Y1, ..., Y4 that contain functions that generate the
% regressors Y1 through Y4 from DOI:10.1002/acs.2945 Section 4
syms p1 p2 p3 fd1 fd2 fs1 fs2 alpha
syms q [2 1]
syms q_dot [2 1]
syms q_del [2 1]
syms q_dot_del [2 1]
syms q_dot_dot [2 1]
syms qd [2 1]
syms qd_dot [2 1]
syms qd_dot_dot [2 1]
e = qd - q;
e_dot = qd - q_dot;
theta = [p1 p2 p3 fd1 fd2 fs1 fs2].';
M = @(q) [p1+2*p3*cos(q(2))   , p2+p3*cos(q(2));
          p2+p3*cos(q(2))     , p2            ];
M_dot = diff(M(q),q1)*q_dot(1) + diff(M(q),q2)*q_dot(2);
Vm = @(q,q_dot) [-p3*sin(q(2))*q_dot(2), -p3*sin(q(2))*(q_dot(1)+q_dot(2));
                 p3*sin(q(2))*q_dot(1) , 0                               ];
Fd = diag([fd1,fd2]);
Fs = @(q_dot) [fs1*tanh(q_dot(1));fs2*tanh(q_dot(2))];
Y1theta = M(q)*q_dot_dot + Vm(q,q_dot)*q_dot + Fd*q_dot + Fs(q_dot);
Y2theta = M(q)*(qd_dot_dot + alpha*e_dot) ...
    + Vm(q,q_dot)*(qd_dot + alpha*e) + Fd*q_dot + Fs(q_dot);
Y3theta = M(q)*q_dot - M(q_del)*q_dot_del;
Y4theta = -M_dot*q_dot + Vm(q,q_dot)*q_dot + Fd*q_dot + Fs(q_dot);

Y1 = sym(zeros(P.n,P.L));
Y2 = sym(zeros(P.n,P.L));
Y3 = sym(zeros(P.n,P.L));
Y4 = sym(zeros(P.n,P.L));
for i = 1:P.n
    [C,T]=coeffs(Y1theta(i),[p1 p2 p3 fd1 fd2 fs1 fs2]);
    Y1(i,find(has(theta,T))) = C;
    [C,T]=coeffs(Y2theta(i),[p1 p2 p3 fd1 fd2 fs1 fs2]);
    Y2(i,find(has(theta,T))) = C;
    [C,T]=coeffs(Y3theta(i),[p1 p2 p3 fd1 fd2 fs1 fs2]);
    Y3(i,find(has(theta,T))) = C;
    [C,T]=coeffs(Y4theta(i),[p1 p2 p3 fd1 fd2 fs1 fs2]);
    Y4(i,find(has(theta,T))) = C;
end
matlabFunction(Y1,'Vars',{q,q_dot,q_dot_dot},'File','Y1','Optimize',false);
matlabFunction(Y2,'Vars',{q,q_dot,qd,qd_dot,qd_dot_dot,alpha},'File','Y2','Optimize',false);
matlabFunction(Y3,'Vars',{q,q_dot,q_del,q_dot_del},'File','Y3','Optimize',false);
matlabFunction(Y4,'Vars',{q,q_dot},'File','Y4','Optimize',false);

% Create a matlab file that contains the function that generates the
% desired trajectry and its derivatives
syms t
qd = [(1+10*exp(-2*t)).*sin(t); (1+10*exp(-t)).*cos(3*t)];
matlabFunction(qd,diff(qd,t),diff(diff(qd,t),t),'File','desiredTraj','Optimize',false);