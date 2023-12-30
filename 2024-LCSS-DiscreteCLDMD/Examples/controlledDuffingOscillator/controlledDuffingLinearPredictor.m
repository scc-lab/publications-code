% This is an implementation of the linear predictor from
%
% M. Korda and I. Mezic, “Linear predictors for nonlinear dynamical systems: 
% Koopman operator meets model predictive control,” Automatica, vol. 93, pp. 149–160, 2018
% 
% The linear predictor is used as baseline for comparsions in Experiment 2 of the paper.
%
% © Zachary Morrison

clear all
close all
clc

%% System Dyanamics 
% Duffing Oscillator
alpha = 1;
beta = -1;
delta = 0;
f = @(t,x,u) [x(2,:) ; -delta.*x(2,:)-beta.*x(1,:)-alpha.*x(1,:).^3 + (2 + sin(x(1,:))).*u];
% g = @(x) [0 ; 2 + sin(x(1))];
n = 2; % States
m = 1; % control input 

%% Runge-Kutta Discretization
h = 0.01; % step size
k1 = @(t,x,u) ( f(t, x, u) );
k2 = @(t,x,u) ( f(t + h/2, x + h/2 * k1(t,x,u), u) );
k3 = @(t,x,u) ( f(t + h/2, x + h/2 * k2(t,x,u), u) );
k4 = @(t,x,u) ( f(t + h, x + h * k3(t,x,u), u) ); % They put k1 here
fd = @(t,x,u) ( x + h/6 * (k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)) );

%% Simulations
Nsim = 1;
u = 2*rand(Nsim,1000)-1; % Matrix of random control inputs in [-1, 1].
xo = 2*rand(2,1000)-1; % Random intiial conditions inputs in [-1, 1]^2.
% Data Matrices & Control Matrix
X = [];
Y = [];
U = [];

for i=1:Nsim
    xn = fd(0, xo, u(i,:));
    X = [X xo];
    Y = [Y xn];
    U = [U u(i,:)];
    xo = xn;
end

% Now we need to perform the lifting of the data matrices X and Y
% Xlift contains X1,X2, then 100 thin plate spline radial basis functions 
% Create the splines
psi_x = zeros(100, length(X));
psi_y = zeros(100, length(Y));
z = zeros(102,1);
centvector = zeros(100,2);
z(1,1) = 2;
z(2,1) = -2;
zo = [z(1,1),z(2,1)];

for k = 1:100
    cent = [unifrnd(-1, 1),unifrnd(-1, 1)];
    centvector(k,1:2) = cent; % keep track of the lifting function centers
    for j=1:length(X)
        x = [X(1,j), X(2,j)];
        y = [Y(1,j), Y(2,j)];
        psi_x(k, j) = (norm(x - cent))^2 * log(norm(x - cent));
        psi_y(k, j) = (norm(y - cent))^2 * log(norm(y - cent));
    end
    z(k+2,1) = norm((zo - cent))^2 * log(norm(zo - cent));
end

% %% New Lifted State
Xlift = [X; psi_x];
Ylift = [Y; psi_y];
 
% Find matrices A, B
% Initialize A in R(n x n), B in R(N x m)
% Here n = 200,000 and so m = 1
n = length(Xlift); % Size of state output
m = length(U); % Size of control sequence

% analytical solution
G = Ylift * pinv([Xlift; U]);

% [A, B] = G
A = G(:, 1:length(G)-1);
B = G(:, length(G));

% This gives us both A and B matrices for our linear predictor
 
% Finding C, C = X * inv(Xlift)
C = X * pinv(Xlift);

%% Test code for the A and B matrices of the linear predictor
SimT = 1; %Simulation Time
Nsim = SimT/h;
ts = h;
t = linspace(0,SimT,Nsim);
% Square wave control input with unit magnitude and freq of 3.33 Hz
% usq =@(i)((-1).^(round(i/30)));
ax = [2, -2]'; % Initial Condition
mu = @(x) -2*x(1,:) - 2*x(2,:);
xhat = zeros(2,Nsim);

for i=1:Nsim
    if i == 1
        xhat = ax;
    else
        z = [z, A*z(:,end) + B*mu(xhat(:,i-1))];
        xhat = [xhat, C*z(:,i)];
    end
end

xhat = C*z;
f = @(x) [x(2) ; -delta*x(2)-beta*x(1)-alpha*x(1)^3];
g = @(x) [0 ; 2 + sin(x(1))];
[~,ax] = ode45(@(t,x) f(x)+g(x)*mu(x), t, ax);
ax = ax.';

%% DCLDMD
% Kernels
kT = 100;
k = 100;
l = 1e-6;

K=KernelvvRKHS('Gaussian',k*ones(m+1,1));
KT=KernelRKHS('Gaussian',kT);
[~,~,~,~,dr,fHat] = controlKoopmanDMD(KT,K,X,U,Y,ts,mu,l);

%% Indirect discrete reconstruction
t_pred = 0:ts:SimT;
y = zeros(2,numel(t_pred));
y(:,1) = [2;-2];
y_pred = zeros(2,numel(t_pred));
y_pred(:,1) = [2;-2];
for i=1:numel(t_pred)-1
y_pred(:,i+1) = dr(1,y_pred(:,i));
[~,temp] = ode45(@(t,x) f(x) + g(x) * mu(x),[0,ts],y(:,i));
y(:,i+1) = temp(end,:).';
end

% %Plots
plot(t,ax,'linewidth',2)
hold on
set(gca,'ColorOrderIndex',1)
plot(t,y_pred(:,1:end-1), t, xhat,'--','linewidth',2)
hold off
xlabel('Time (s)')
set(gca,'fontsize',16)
legend('$x_{1}(t)$','$x_{2}(t)$','$\hat{x}_1(t)$','$\hat{x}_{2}(t)$','$x_{p,1}(t)$','$x_{p,2}(t)$',...
'interpreter','latex','fontsize',16,'location','southeast')
% DiscreteReconCompare = [t_pred(:,1:end-1)',y_pred(:,1:end-1)', xhat', ax']; % data for tikzplot
% save('DiscreteReconCompare.dat','DiscreteReconCompare','-ascii');

%%Plotting the simulated Solution 
% figure;
% subplot(2,2,1)
% plot([1:Nsim]*h, ax(1,:), 'b')
% ylim([-1.25 1.25])
% xlabel('time [s]','fontsize',12,'interpreter','latex')
% ylabel('$x_1$','fontsize',14,'interpreter','latex')
% 
% subplot(2,2,2)
% plot([1:Nsim]*h, ax(2,:), 'r')
% ylim([-1.25 1.25])
% xlabel('time [s]','fontsize',12,'interpreter','latex')
% ylabel('$x_2$','fontsize',14,'interpreter','latex')


%% Plotting the Predictor Solution
% subplot(2,2,3)
% plot([1:Nsim]*h, xhat(1,:), '--b')
% ylim([-1.25 1.25])
% xlabel('time [s]','fontsize',12,'interpreter','latex')
% ylabel('$\hat{x}_{1}$','fontsize',14,'interpreter','latex')
% 
% subplot(2,2,4)
% plot([1:Nsim]*h, xhat(2,:), '--r')
% ylim([-1.25 1.25])
% xlabel('time [s]','fontsize',12,'interpreter','latex')
% ylabel('$\hat{x}_{2}$','fontsize',14,'interpreter','latex')

%% Together
lw = 2;

% hold on
% figure
% plot([1:Nsim]*h, ax(1,:), [1:Nsim]*h, xhat(1,:), '--b','linewidth',lw)
% axis([1 SimT min(xhat(1,:))-0.15 max(xhat(1,:))+0.15])
% title('Predictor comparison - $x_1$','interpreter','latex'); 
% xlabel('Time [s]','interpreter','latex');
% set(gca,'fontsize',20)
% LEG = legend('True','Koopman');
% set(LEG,'interpreter','latex')
% 
% figure
% plot([1:Nsim]*h, ax(2,:), [1:Nsim]*h, xhat(2,:), '--r','linewidth',lw)
% axis([1 SimT min(xhat(2,:))-0.15 max(xhat(2,:))+0.15])
% title('Predictor comparison - $x_2$','interpreter','latex'); 
% xlabel('Time [s]','interpreter','latex');
% set(gca,'fontsize',20)
% LEG = legend('True','Koopman');
% set(LEG,'interpreter','latex')
