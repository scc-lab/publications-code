% Problem parameters
P.n=1;
P.m=1;
P.R = 1;
P.W_Q_star = [0;1];
P.Q = diag(P.W_Q_star);
P.W_V_star = [pi/2;1;1];
P.W_star=[P.W_V_star;P.W_Q_star];
P.theta_star = [-1;-5/2;4];
P.p_0 = ones(P.n,1);
P.q_0 = ones(P.n,1);
P.G = 3;

% Parameter Estimator
PE.L = numel(P.theta_star);
PE.theta_hat_0 = zeros(size(P.theta_star));
PE.Gamma_0 = 1*diag([1;100;1]);
PE.T1 = 1;
PE.T2 = 0.8;
PE.N = 150;
PE.k_theta = 2/PE.N;
PE.beta = 0.1;
PE.t_step = 0.01;
PE.tau = 0.1;

% State Estimator
SE.q_hat_0 = zeros(P.n,1);
SE.k = 15;
SE.alpha = 15;
SE.beta = 15;

% IRL
IRL.N = 150;
IRL.alpha = 0.01/IRL.N;
IRL.beta = 0.6;
IRL.W_hat_0 = zeros(numel(P.W_star),1);
IRL.Gamma_W_0 = 0.001*diag([1;1;1;1;0.1]);
IRL.tau = 0.1;
IRL.num_weights = size(P.W_star,1);
IRL.t_step = 0.01;