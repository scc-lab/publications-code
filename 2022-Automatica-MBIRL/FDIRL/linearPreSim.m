% Problem parameters
P.n=1;
P.m=1;
P.A1 = -0.5;
P.A2 = -0.5;
P.A = [zeros(P.n) eye(P.n);P.A1 P.A2];
P.Ad = [zeros(P.n) eye(P.n);-2*eye(P.n) zeros(P.n)];
P.B1 = ones(P.n,P.m);
P.B = [zeros(P.n,P.m);P.B1];
P.B_plus = pinv(P.B);
P.R = 50*eye(P.m);
P.Q = diag([1.1 3]);
[P.kp,P.P_opt,~]  = lqr(P.A,P.B,P.Q,P.R);
P.W_star=[P.P_opt(1,1);P.P_opt(2,2);2*P.P_opt(1,2);P.Q(1,1);P.Q(2,2)];
P.theta_star=[reshape(P.A1,P.n^2,1);reshape(P.A2,P.n^2,1)];
P.p_0 = ones(P.n,1);
P.q_0 = ones(P.n,1);
P.x_d_0 = [0;2];

% Parameter Estimator
PE.L = numel(P.theta_star);
PE.theta_hat_0 = zeros(size(P.theta_star));
PE.Gamma_0 = eye(numel(P.theta_star));
PE.T1 = 0.4;
PE.T2 = 0.2;
PE.Lags = [PE.T1;PE.T2;PE.T1+PE.T2];
PE.N = 50;
PE.k_theta = 0.5/PE.N;
PE.beta = 2;
PE.t_step = 0.01;

% State Estimator
SE.q_hat_0 = zeros(P.n,1);
SE.k = 5;
SE.alpha = 5;
SE.beta = 5;

% Feedback Estimator
FE.K = 2;
FE.M = 40;
FE.beta_u = 1;
FE.alpha_u = 1;
FE.tau = 0.1;
FE.Gamma_Wu_0 = 0.0002*eye(2);
FE.Wu_0 = 0*ones(2,1);
FE.t_step = 0.01;

% IRL
IRL.feedback_driven = 1;
IRL.N = 10;
IRL.alpha = 0.1/IRL.N;
IRL.beta = 0.2;
IRL.W_hat_0 = zeros(numel(P.W_star),1);
IRL.Gamma_W_0 = 1*eye(numel(P.W_star));
IRL.tau = 1;
IRL.max_e = 5;
IRL.num_weights = size(P.W_star,1);
IRL.t_step = 0.01;