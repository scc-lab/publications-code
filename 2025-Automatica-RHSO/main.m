% Code Saved 12/10/2025
% Created by Jared Town and Rushikesh Kamalapurkar

rng(1) % This reproduces results in the paper, comment it out to randomize
global dellta
global glob_sig_hat
global glob_sig_u
global X_store_rank
global t_batch
global Z_store_rank
global sig_hat_glob
global sig_u_glob 

% Sim Steps
P.t_step = 0.01;
P.t_end = 300;
P.eps = 0.001; % regularization parameter for the RHSO
P.t_dwell = 0.05; % dwell time for history stack data recording
P.t_purge = 2; % dwell time for purging the history stack
P.t_rank = 0.5; % time period to check FI conditions

% System model
unique = 0; % Set this to 1 for an IRL problem with a unique solution
% Toy system
A_diag = diag([1 2 3]);
B_diag = diag([1 1 1]);
C_diag = eye(3);
% Transformation matrix
T = [1 2 -1; 
     -1 3 4;
     1 2 -3];
% Transformed system
P.A = T*A_diag*inv(T);
P.B = T*B_diag;
P.C = C_diag*inv(T);
P.n = size(P.A, 1); % Number of states
P.m = size(P.B, 2); % Number of inputs
P.x_0 = ones(P.n,1); % Initial condition
% Adding this value makes the IRL problem admit a unique solution
if unique == 0
    P.A(1,1)=1;
end
% Number of outputs
P.L = size(P.C, 1);

%%% Cost function matrices
Q = diag([7, 4, 8]);
P.Q = inv(T).'*Q*inv(T);
P.R = diag([1, 4, 7]);

P.M = P.m;%(P.m+1)/2; % Numer of unknowns in R (Assumed diagonal)
P.P_v = P.n*(P.n+1)/2; % Number of unknowns in S
P.P_q = P.n*(P.n+1)/2; % Number of unknowns in Q
P.N = (P.P_v + P.P_q + P.M-1); % Total number of unknowns

P.r1 = P.R(1, 1); % First control penalty assumed to be known
[P.k_p,P.Opt_S] = lqr(P.A,P.B,P.Q,P.R); % Optimal controller

% Ideal Weights
P.W_V_Star = zeros(P.P_v, 1);
k = 1;
for i = 1:P.n
    for j = (1+i-1):(P.n)
        P.W_V_Star(k) = P.Opt_S(i,j);
        k = k+1;
    end
end
P.W_Q_Star = zeros(P.P_q, 1);
k = 1;
for i = 1:P.n
    for j = (1+i-1):(P.n)
        P.W_Q_Star(k) = P.Q(i,j);
        k = k+1;
    end
end
P.W_R_Star = zeros(P.m-1,1);
for i = 1:P.m-1
    P.W_R_Star(i) = P.R(i+1,i+1);
end
P.W_Star = [P.W_V_Star;P.W_Q_Star;P.W_R_Star];

% Estimation gains
P.K_3 = place(P.A',P.C',[-0.1 -1.5 -2]).';
P.K_4 = eye(P.P_v+P.P_q+P.M-1);

% Initial conditions
% P.W_hat_0 = 3*ones(P.P_v+P.P_q+P.M-1,1); % Constant initialization
P.W_hat_0 = randn(P.P_v+P.P_q+P.M-1,1); % Random initialization
P.x_hat_0 = 2*ones(P.n,1);

% Memory allocation for basis functions
P.sigma_V = zeros(P.P_v, 1);
P.sigma_R1 = zeros(1, P.M-1);
P.sigma_Q = zeros(P.P_q, 1);

% Frequencies and phases for excitation signal
P.freq_mag = 1;
P.freq = [0.001*rand(5,P.m); 0.01*rand(5,P.m); 0.1*rand(5,P.m); 1*rand(5,P.m); 10*rand(5,P.m)];
P.phase = pi*rand(5*5,P.m);

%%%% Running the main simulation %%%%
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,z]=ode45(@(t,z) closedLoop(t,z,P),0:P.t_step:P.t_end,[P.x_0;P.x_hat_0;P.W_hat_0], opts);

%%% Pulling out final values to do analysis with
W_hat_ALL = z(:, P.n*2+1:end);
W_hat_V_ALL = W_hat_ALL(:, 1:P.P_v);
W_hat_Q_ALL = W_hat_ALL(:, P.P_v+1:P.P_v+P.P_q);
W_hat_R_prime_ALL = W_hat_ALL(:, P.P_v+P.P_q+1:end);
W_hat_R_ALL = [repmat(P.r1, size(W_hat_R_prime_ALL,1), 1), W_hat_R_prime_ALL];
x_ALL = z(:, 1:P.n);
x_hat_ALL = z(:, P.n+1:P.n*2);
W_hat = z(end,P.n*2+1:end);
W_hat_V = W_hat(1:P.P_v);
W_hat_Q = W_hat(P.P_v+1:P.P_v+P.P_q);
W_hat_R = [P.r1, W_hat(P.P_v+P.P_q+1:end)];

R_hat = diag(W_hat_R);
Q_t = zeros(P.n);
k=1;
for i = 1:P.n
    for j = (1+i-1):(P.n)
        Q_t(i,j) = W_hat_Q(k);
        k = k+1;
    end
end
Q_hat=Q_t'+tril(Q_t',-1)';
S_hat = zeros(P.n);
k=1;
for i = 1:P.n
    for j = (1+i-1):(P.n)
        S_hat(i,j) = W_hat_V(k);
        k = k+1;
    end
end
S_hat=S_hat'+tril(S_hat',-1)';
k_p_hat = R_hat\P.B.'*S_hat;
% Check if Q_hat, R_hat, and S_hat are compatible
[X, L, G] = care(P.A, P.B, Q_hat, R_hat);
diff_S = S_hat-X;


k_p_hat_ALL = zeros(size(W_hat_ALL, 1), numel(P.k_p));
HSO_ARE_ALL = zeros(size(W_hat_ALL, 1), P.n^2);
HSO_HJB_ALL = zeros(size(W_hat_ALL, 1), 1);
HJB_ALL = zeros(size(W_hat_ALL, 1), 1);
for n = 1:size(W_hat_ALL, 1)
    % Generating Q, R, and S matrices from their respective basis vectors 
    % at each time instance
    W_V_t = W_hat_V_ALL(n, :);
    W_Q_t = W_hat_Q_ALL(n, :);
    x_hat_t = x_hat_ALL(n, :).';
    x_t = x_ALL(n, :).';
    
    R_hat_t = diag(W_hat_R_ALL(n, :));
    S_t = zeros(P.n);
    Q_t = zeros(P.n);
    k=1;
    % Constructing the symmetric S and Q matrices.
    for i = 1:P.n
        for j = (1+i-1):(P.n)
            S_t(i,j) = W_V_t(k);
            Q_t(i,j) = W_Q_t(k);
            k = k+1;
        end
    end
    S_t=S_t'+tril(S_t',-1)';
    Q_t=Q_t'+tril(Q_t',-1)';
    
    k_p_hat_t = R_hat_t\P.B.'*S_t;
    k_p_hat_ALL(n, :) = k_p_hat_t(:);
    %HSO_ARE_temp = P.A.'*HSO_S_temp + HSO_S_temp*P.A*P.B*inv(HSO_R_temp)*P.B.'*HSO_S_temp+HSO_Q_temp;
    % I am re-writing the above equation for easier readability
    Q = Q_t;
    R = R_hat_t;
    S = S_t;
    A = P.A;
    B = P.B;
    HSO_ARE_t = A.'*S + S*A - S*B*inv(R)*B.'*S + Q;
    HSO_HJB_ALL(n) = x_hat_t.'*(A.'*S + S*A - S*B*inv(R)*B.'*S + Q)*x_hat_t;
    HJB_ALL(n) = x_t.'*(A.'*P.Opt_S + P.Opt_S*A - P.Opt_S*B*inv(P.R)*B.'*P.Opt_S + P.Q)*x_t;
    HSO_ARE_ALL(n,:) = HSO_ARE_t(:);
end

%%% Running the learned system %%%%
[k_p_hat_lqr,S_hat_lqr] = lqr(P.A,P.B,Q_hat,R_hat);
newClosedLoop = @(t,x,P) xdot(t,x,-k_p_hat_lqr*x+uexc(t,P),P);
[~,new_x]=ode45(@(t,x) newClosedLoop(t,x,P),0:P.t_step:P.t_end,P.x_0, opts);

%%% Checking FI conditions
stem_value=[];
temp = [];
for i = 1:14:size(glob_sig_hat,2)-1
    space = glob_sig_hat(:,i:i+14-1);
    if space == 1
        continue
    end
    space_rank = rank(space);
    space_new_rank = rank([space, glob_sig_u(:,(i+14-1)/14)]);
    if space_rank == space_new_rank
        stem_value = [stem_value; 1];
    else
        temp = [temp; (i+14-1)/14];
        stem_value = [stem_value; 0];
    end
end

x_stem_values = [];
for i = 1:size(X_store_rank, 1)
    if X_store_rank(i) == P.n
        x_stem_values = [x_stem_values; 1];
    else
       x_stem_values = [x_stem_values; 0];
    end
end

z_stem_values = [];
for i = 1:size(Z_store_rank)
    if Z_store_rank(i) == P.n*(P.n+1)/2
        z_stem_values = [z_stem_values; 1];
    else
        z_stem_values = [z_stem_values; 0];
    end
end

% Ridge regression test
Sig_u = sig_u_glob;
Sig_hat = sig_hat_glob;
W_hat_ridge = (Sig_hat.'*Sig_hat+P.eps*eye(size(Sig_hat.'*Sig_hat)))\(Sig_hat.'*Sig_u);
W_hat_V_ridge = W_hat_ridge(1:P.P_v);
W_hat_R_ridge = [P.r1, W_hat_ridge(P.P_v+P.P_q+1:end).'];
R_hat_ridge = diag(W_hat_R_ridge);
W_hat_Q_ridge = W_hat_ridge(P.P_v+1:P.P_v+P.P_q);
Q_t = zeros(P.n);
k=1;
for i = 1:P.n
    for j = (1+i-1):(P.n)
        Q_t(i,j) = W_hat_Q_ridge(k);
        k = k+1;
    end
end
Q_hat_ridge=Q_t'+tril(Q_t',-1)';

[k_p_hat_ridge_lqr,S_hat_ridge_lqr] = lqr(P.A,P.B,Q_hat_ridge,R_hat_ridge);
newClosedLoop_ridge = @(t,x,S) xdot(t,x,-k_p_hat_ridge_lqr*x+uexc(t,S),S);
[~,new_x_ridge]=ode45(@(t,x) newClosedLoop_ridge(t,x,P),0:P.t_step:P.t_end,P.x_0, opts);

%%%% Data for plots in the paper
if unique == 1
    uniqueString = '_unique';
else
    uniqueString = '';
end
downsam = 200;
temp = [t, vecnorm((new_x-z(:,1:P.n)).').'];
temp = temp(1:downsam:size(temp,1),:);
save(['x_norm_error' uniqueString '.dat'],'temp','-ascii');

temp = [t,vecnorm(z(:,P.n*2+1+P.P_v:P.n*2+P.P_v+P.P_q).'-P.W_Q_Star).'];
temp = temp(1:downsam:size(temp,1),:);
save(['Q_norm_error' uniqueString '.dat'],'temp','-ascii');

temp = [t, vecnorm(z(:,end-(P.m-2):end).'-P.W_R_Star,2,1).'];
temp = temp(1:downsam:size(temp,1),:);
save(['R_norm_error' uniqueString '.dat'],'temp','-ascii');

temp = [t,vecnorm((k_p_hat_ALL - P.k_p(:).').').'];
temp = temp(1:downsam:size(temp,1),:);
save(['kp_error' uniqueString '.dat'],'temp','-ascii');

downsam_batch = 3;
temp = [t_batch.', vecnorm(dellta).'];
temp = temp(1:downsam_batch:size(temp,1),:);
save(['Delta_norm' uniqueString '.dat'],'temp','-ascii');

temp = [t_batch.', stem_value];
temp = temp(1:downsam_batch:size(temp,1),:);
save(['Sigma_u_condition' uniqueString '.dat'],'temp','-ascii');

temp = [t_batch.', x_stem_values];
temp = temp(1:downsam_batch:size(temp,1),:);
save(['x_span_condition' uniqueString '.dat'],'temp','-ascii');

temp = [t_batch.', z_stem_values];
temp = temp(1:downsam_batch:size(temp,1),:);
save(['symmetric_span_condition' uniqueString '.dat'],'temp','-ascii');

temp = [t, vecnorm((new_x_ridge - z(:,1:P.n)).').'];
temp = temp(1:downsam:size(temp,1),:);
save(['x_norm_error_ridge' uniqueString '.dat'],'temp','-ascii');


% Test: use all data to compute k
% x_ALL = x_ALL.'; %Transposing this to make the math easier on the eyes
% x_hat_ALL = x_hat_ALL.';
% u_ALL = -P.k_p*x_ALL;
% k_check = (-u_ALL*x_ALL.')/(x_ALL*x_ALL.');
% P.k_p-k_check;
% k_check_hat = (-u_ALL*x_hat_ALL.')/(x_hat_ALL*x_hat_ALL.');
% P.k_p-k_check_hat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plotting Original Simulation measurements
% % States
% figure
% plot(t,z(:,1:P.n)); title('Measured trajectory')
% 
% % State estimation error
% figure
% semilogy(t,abs(z(:,1:P.n)-z(:,P.n+1:P.n*2))); title('State estimation error')
% 
% figure
% plot(t,new_x); title('Optimal trajectory for learned cost function')
% figure
% plot(t,new_x - z(:,1:P.n)); title('Trajectory prediction error')
% xlabel('time (s)')
% figure
% semilogy(t, vecnorm(new_x.'-z(:,1:P.n).'), LineWidth=1.5) % Trajectory Difference
% ylim([0 1])
% xlabel('Time (s)')
% ylabel('$||(X_{Expert}-X_{Learner})^{T}||$', Interpreter='latex', FontSize=15)
% 
% figure
% semilogy(t_batch, vecnorm(dellta),LineWidth=1.5); % Delta Value
% ylabel('$||\Delta||$', Interpreter='latex', FontSize=15)
% xlabel('Time (s)', FontSize=12)
% 
% figure
% plot(t,z(:,P.n*2+1+P.P_v:P.n*2+P.P_v+P.P_q)-P.W_Q_Star.',LineWidth=1.5); % Q Difference
% xlabel('Time (s)', FontSize=12)
% ylabel('$Q_{Expert}-Q_{Learner}$', Interpreter='latex', FontSize=15)
% 
% figure
% plot(t,z(:,end-(P.m-2):end)-P.W_R_Star',LineWidth=1.5); % R Difference
% xlabel('Time (s)', FontSize=12)
% ylabel('$R_{Expert}-R_{Learner}$', Interpreter='latex', FontSize=15)
% 
% figure
% plot(t,z(:,P.n*2+1+P.P_v:P.n*2+P.P_v+P.P_q)-P.W_Q_Star.', 'r',LineWidth=1.5); % Q Difference
% hold on
% plot(t,z(:,end-(P.m-2):end)-P.W_R_Star', 'b',LineWidth=1.5); % R Difference
% hold off
% xlabel('Time (s)', FontSize=12)
% ylabel('$Q,R_{Expert}-Q,R_{Learner}$', Interpreter='latex', FontSize=15)
% 
% figure
% plot(t,vecnorm(z(:,P.n*2+1+P.P_v:P.n*2+P.P_v+P.P_q).')-vecnorm(P.W_Q_Star.'), 'r',LineWidth=1.5); % Q Difference
% hold on
% % Taking the vector norm of a vector yeilds a scalar, conditional is to fix that.
% if size(P.W_R_Star, 1) == 1
%     plot(t,z(:,end-(P.m-2):end).'-P.W_R_Star', 'b',LineWidth=1.5); % R Difference
% else
%     plot(t,vecnorm(z(:,end-(P.m-2):end).'-P.W_R_Star), 'b',LineWidth=1.5); % R Difference
% end
% hold off
% xlabel('Time (s)', FontSize=12)
% ylabel('$||Q,R_{Expert}||-||Q,R_{Learner}||$', Interpreter='latex', FontSize=15)
% legend('$||Q_{Diff}||$', '$||R_{Diff}||$', Interpreter='latex')
% 
% figure
% semilogy(t,vecnorm((k_p_hat_ALL - P.k_p(:).').'),LineWidth=1.5) % Gain Difference
% xlabel('Time (s)', FontSize=12)
% ylabel('$||Kp_{Expert}-Kp_{Learner}||$', Interpreter='latex', FontSize=15)
% 
% delta = dellta(:, weird);
% time_delta = t_batch(:,weird);
% 
% stem(t_batch, stem_value)
% ylabel('$\Sigma_u=$Range$(\hat{\Sigma})$', Interpreter='latex', FontSize=15)
% xlabel('Time (s)', FontSize=12)
% figure
% plot(t,new_x); title('Optimal trajectory for learned cost function - Ridge Regression')
% figure
% plot(t,vecnorm((new_x_ridge - z(:,1:P.n)).')); title('Trajectory prediction error - Ridge Regression')
% xlabel('time (s)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_exc = uexc(t,P)
    u_exc = zeros(P.m,1);
    for i=1:size(P.freq,1)
        u_exc = u_exc + P.freq_mag*sin(2*pi*P.freq(i,:)*t.' + P.phase(i,:)).';
    end
end

function xD = xdot(t,x,u,P)
    xD = P.A*x + P.B*u;
end

function [out1, out2, out3] = history_deltas(deltas, u, t, P, x_hat)
% History stack management
    persistent time_last_stepped
    persistent s_hat_current
    persistent s_u_current
    persistent s_hat
    persistent s_u
    persistent jdx
    persistent t_last_update
    persistent X_store
    persistent x_idx
    global flag
    global X_store_current
    if t == 0
        time_last_stepped = -5;
        X_store_current = zeros(P.n,P.N);
        X_store = zeros(P.n, P.N);
        x_idx = 1;
        s_hat_current = zeros(P.N+P.N*P.m,P.P_v+P.P_q+P.M-1);
        s_u_current = zeros(P.N+P.N*P.m,1);
        s_hat = zeros(P.N+P.N*P.m,P.P_v+P.P_q+P.M-1);
        s_u = zeros(P.N+P.N*P.m,1);
        jdx = 1;
        t_last_update = t;
        flag = false; % Flag is used to populate the first history stack without checking rank
    end
    new_cond = 1e60;
    % Increments in steps of 2 rows
    if jdx < P.N+P.N*P.m && (t-time_last_stepped) >= P.t_dwell
        % Storing the first full set of data
        s_hat(jdx:jdx+P.m,:) = deltas;
        s_u(jdx:jdx+P.m,:) = u;
        X_store(:,x_idx) = x_hat;
        jdx = jdx + (P.m+1);
        x_idx = x_idx+1;
        time_last_stepped = t;
        if flag == false && jdx == (P.N+P.N*P.m+1)
        % populate the first history stack without checking conditioning
            s_hat_current = s_hat;
            s_u_current = s_u;
            X_store_current = X_store;
            flag = true; 
        end
    elseif  jdx == P.N+P.N*P.m+1 && t-time_last_stepped >= P.t_dwell
        %%% condition number minimization to populate the temp history stack
        time_last_stepped = t;
        REG = P.eps*eye(size(s_hat.'*s_hat));
        min_store = zeros(P.N,1);
        current_sum = s_hat.'*s_hat + REG;
        current_min = cond(current_sum);
        for idx=1:(P.m+1):(P.N+P.N*P.m) % Number of values saved in the history stack
           temp_sum = current_sum - s_hat(idx:idx+P.m,:).'*s_hat(idx:idx+P.m,:) + deltas.'*deltas;
           temp_min = cond(temp_sum);
           v = round(1/(P.m+1)*idx + (1-1/(P.m+1)));
           min_store(v) = temp_min;
        end
        [new_min, index] = min(min_store);
        if new_min < current_min
           s_hat(index*(P.m+1)-P.m:(P.m+1)*index,:) = deltas;
           new_sum = s_hat.'*s_hat + REG;
           new_cond = cond(new_sum);
           s_u(index*(P.m+1)-P.m:(P.m+1)*index,:) = u;
           X_store(:,index) = x_hat;
        end
    end
    %%% Update in-use history stack and purge the temp history stack
    if jdx == P.N+P.N*P.m+1 && t-t_last_update >  P.t_purge || new_cond < 1e+5
        X_store_current = X_store;
        s_hat_current = s_hat;
        s_u_current = s_u;
        jdx = 1;
        x_idx = 1;
        s_hat = zeros(P.N+P.N*P.m,P.P_v+P.P_q+P.M-1);
        s_u = zeros(P.N+P.N*P.m, 1);
        X_store = zeros(P.n, P.N);
        t_last_update = t;
    end
    out1 = s_hat_current;
    out2 = s_u_current;
end

% Linear System
function zDot = closedLoop(t, Z, P)
    global sig_hat_glob
    global sig_u_glob
    global flag
    global dellta
    global glob_sig_hat
    global glob_sig_u
    persistent last_t 
    global X_store_current
    global X_store_rank
    %global BatchWhat
    global t_batch
    global Z_store_rank

    x = Z(1:P.n);
    x_hat = Z(P.n+1:P.n*2);
    W_hat = Z(P.n*2+1:end);
    u = -P.k_p*x;
    % excitation signal
    u_exc = zeros(P.m,1);
    for i=1:size(P.freq,1)
        u_exc = u_exc + P.freq_mag*sin(2*pi*P.freq(i,:)*t.' + P.phase(i,:)).';
    end

    x_dot = P.A*x + P.B*(u+u_exc);
    y_prime = P.C*x;
    
    % Generating sigma_V
    P.sigma_V = Sigma_V(x_hat, P);
    
    % Generating sigma_Q
    P.sigma_Q = Sigma_Q(x_hat, P);
    
    % Generating sigma_R1
    P.sigma_R1_ = u.^2.';
    P.sigma_R1_(:,1) = [];
    
    % Generating sigma_R2
    sigma_R2_ = Sigma_R_Diag(u, P);

    % Getting Values for Gradient of Sigma V
    if ~isfile("Grad_Sigma_V.m")
        Grad_Sigma_V_gen(P);
    end
    grad_sigma_V = Grad_Sigma_V(x);
    
    sigma_delta = [(P.A*x_hat + P.B*u).'*grad_sigma_V.', P.sigma_Q.', P.sigma_R1_];
    sigma_Delta_u = [P.B.'*grad_sigma_V.', zeros(P.m,P.P_q), 2*sigma_R2_];
    sigma_u = [-u(1)^2*P.r1; -2*u(1)*P.r1; zeros(P.m-1,1)];
    [sig_hat, sig_u] = history_deltas([sigma_delta; sigma_Delta_u], sigma_u, t, P, x_hat);
    if flag
        K_4 = P.K_4/(sig_hat.'*sig_hat+eye(size(sig_hat.'*sig_hat))*P.eps)*sig_hat.';
    else
        K_4 = zeros(P.P_v+P.P_q+P.M-1,P.N+P.N*P.m);
    end
    K = [P.K_3                         zeros(P.n,P.N+P.N*P.m);
         zeros(P.P_v+P.P_q+P.M-1,P.L)  K_4];
    innovation = [y_prime; sig_u] - [P.C*x_hat; sig_hat*W_hat];
    x_W_hat_dot = [P.A*x_hat + P.B*(u+u_exc); zeros(P.P_q+P.P_v+P.M-1,1)] + K * innovation;
    zDot = [x_dot;x_W_hat_dot];

    % Code for data storage to check the FI conditions
    if t == 0
        dellta = [];
    end
    if flag && isempty(dellta)
        dellta = sig_u-sig_hat*W_hat;
        glob_sig_hat = sig_hat;
        glob_sig_u = sig_u;
        last_t = t;
        t_batch = t;
        Z = [];
        for i = 1:P.N
            z = [];
            for j = 1:P.n
                for k = j:P.n
                    a = X_store_current(j,i);
                    b = X_store_current(k,i);
                z = [z; a*b];
                end
            end
            Z = [Z,z];
        end
        Z_store_rank = rank(Z);
        X_store_rank = rank(X_store_current);
    elseif flag && t-last_t>P.t_rank
        t_batch = [t_batch t];
        dellta = [dellta sig_u-sig_hat*W_hat];
        glob_sig_hat = [glob_sig_hat sig_hat];
        glob_sig_u = [glob_sig_u sig_u];
        last_t = t;
        Z = [];
        for i = 1:P.N
            z = [];
            for j = 1:P.n
                for k = j:P.n
                    a = X_store_current(j,i);
                    b = X_store_current(k,i);
                z = [z; a*b];
                end
            end
            Z = [Z,z];
        end
        Z_store_rank = [Z_store_rank; rank(Z)];
        X_store_rank = [X_store_rank; rank(X_store_current)];
    end

    sig_hat_glob = sig_hat;
    sig_u_glob = sig_u;
end

% Basis functions
function sigma_Q = Sigma_Q(x,P)
  sigma_Q = zeros(P.P_q,1);
  k = 1;
  for i = 1:P.n
    for j = i:P.n
      if i == j
        sigma_Q(k) = x(i)^2;
        k = k+1;
      else
        sigma_Q(k) = 2*x(i)*x(j);
        k = k+1;
      end
    end
  end
end

function sigma_Q = Sigma_Q_Reduced_Diag(x)
  sigma_Q = [x(1:3);x(9)].^2;
end

function sigma_R = Sigma_R_Diag(u,P)
  sigma_R = zeros(P.m);
  for i = 1:P.m
    sigma_R(i,i) = u(i);
  end
  sigma_R(:,1) = [];
end

function [sigma_V] = Sigma_V(x,P)
  sigma_V = zeros(P.P_v,1);
  k = 1;
  for i = 1:P.n
    for j = i:P.n
      if i == j
        sigma_V(k) = x(i)^2;
        k = k+1;
      else
        sigma_V(k) = 2*x(i)*x(j);
        k = k+1;
      end
    end
  end
end

function Grad_Sigma_V_gen(P)
    k = 1;
    syms x [P.n,1] real
    sigma_V = sym(zeros(P.P_v,1));
    for i = 1:P.n
      for j = i:P.n
        if i == j
          sigma_V(k) = x(i)^2;
          k = k+1;
        else
          sigma_V(k) = 2*x(i)*x(j);
          k = k+1;
        end
      end
    end
    Grad_Sigma_V = jacobian(sigma_V,x);
    matlabFunction(Grad_Sigma_V,File="Grad_Sigma_V",Vars={x});
end
