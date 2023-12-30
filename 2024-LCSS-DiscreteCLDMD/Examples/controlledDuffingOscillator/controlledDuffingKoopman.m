% This script generates trajectories of a controlled Duffing oscillator
% uses control Liouville DMD to generate a predictive model for a given 
% feedback controller.
%
% © Rushikesh Kamalapurkar and Zachary Morrison
function controlledDuffingKoopman()

addpath('../../lib');

%% Generate Trajectories
n = 2; % Number of dimensions that f maps from/to
m = 1; % Dimensions of the controller
alpha = 1;
beta = -1;
delta = 0;
f = @(x) [x(2) ; -delta*x(2)-beta*x(1)-alpha*x(1)^3];
g = @(x) [0 ; 2 + sin(x(1))];
IV_selection = 'grid'; 
samp_min = -3;
samp_max = 3;
if strcmp(IV_selection,'grid')
    pointsPerDim = 15;
    XDim = linspace(samp_min,samp_max,pointsPerDim);
    [XX,YY] = meshgrid(XDim,XDim);
    X = [XX(:) YY(:)].';
    M = size(X,2);
elseif strcmp(IV_selection,'random')
    % Get TotalTrajectories random IV's.
    M = 100;
    X = samp_min + (samp_max - samp_min)*rand(n, M);
elseif strcmp(IV_selection,'halton')
    M = 100;
    % Get TotalTrajectories halton sequence
    haltonseq = @(n,d) net(haltonset(d),n);
    halton = haltonseq(M, n);
    X = samp_min + (samp_max - samp_min)*halton.';
else
    error('Unknown IV selection mode %s', IV_selection)
end
ts = 0.01;
U = -2+4*rand(1,M);
Y=zeros(size(X));
for i = 1:M
    F = @(x,u) f(x) + g(x) * u; % The update function
    [~,y] = ode45(@(t,x) F(x,U(:,i)),[0,ts],X(:,i));
    Y(:,i) = y(end,:).';
end
%% Kernels
kT = 10;
k = 10;
l = 1e-6;

K=KernelvvRKHS('Gaussian',k*ones(m+1,1));
KT=KernelRKHS('Gaussian',kT);

%% Feedback controller
mu = @(x) -2*x(1,:) - 2*x(2,:);

%% CLDMD
[~,~,~,~,dr,fHat] = ControlKoopmanDMD(KT,K,X,U,Y,ts,mu,l);

%% Indirect reconstruction
% x0 = [2;-2];
% t_pred = 0:0.1:10;
% [~,y_pred] = ode45(@(t,x) fHat(x),t_pred,x0);
% [~,y] = ode45(@(t,x) f(x) + g(x) * mu(x),t_pred,x0);

% Plots
% plot(t_pred,y,'linewidth',2)
% hold on
% set(gca,'ColorOrderIndex',1)
% plot(t_pred,y_pred,'--','linewidth',2)
% hold off
% xlabel('Time (s)')
% set(gca,'fontsize',16)
% legend('$x_1(t)$','$x_2(t)$','$\hat{x}_1(t)$','$\hat{x}_2(t)$',...
% 'interpreter','latex','fontsize',16,'location','southeast')

%% Indirect discrete reconstruction
t_pred = 0:ts:6;
y = zeros(n,numel(t_pred));
y(:,1) = [2;-2];
y_pred = zeros(n,numel(t_pred));
y_pred(:,1) = [2;-2];
for i=1:numel(t_pred)-1
y_pred(:,i+1) = dr(1,y_pred(:,i));
[~,temp] = ode45(@(t,x) f(x) + g(x) * mu(x),[0,ts],y(:,i));
y(:,i+1) = temp(end,:).';
end
% Plots
% plot(t_pred,y,'linewidth',2)
% hold on
% set(gca,'ColorOrderIndex',1)
% plot(t_pred,y_pred,'--','linewidth',2)
% hold off
% xlabel('Time (s)')
% set(gca,'fontsize',16)
% legend('$x_1(t)$','$x_2(t)$','$\hat{x}_1(t)$','$\hat{x}_2(t)$',...
% 'interpreter','latex','fontsize',16,'location','southeast')
% DiscreteReconCompare = [t_pred',y_pred',y']; % data for tikzplot
% DiscreteReconError = [t_pred',y_pred'-y']; 
% save('DiscreteReconCompare.dat','DiscreteReconCompare','-ascii');
% save('DiscreteReconError.dat','DiscreteReconError','-ascii');

%% Vector Field Plot
XDimeval = linspace(-2,2,9);
[XX, YY] = meshgrid(XDimeval,XDimeval);
IVeval = [XX(:) YY(:)].';
x_dot_hat_at_x0 = [];
x_dot_at_x0 = [];
for i=1:size(IVeval,2)
    x0=IVeval(:,i);
    x_dot_hat_at_x0 = [x_dot_hat_at_x0, fHat(x0)];
    x_dot_at_x0 = [x_dot_at_x0, f(x0)+g(x0)*mu(x0)];
end
%max(max(abs(x_dot_at_x0 - x_dot_hat_at_x0)))
figure
subplot(2,2,1);
surf(XX,YY,reshape(x_dot_hat_at_x0(1,:),9,9))
xlabel('$x_1$','interpreter','latex','fontsize',16)
ylabel('$x_2$','interpreter','latex','fontsize',16)
zlabel('$\left(\hat{f}(x) + \hat{g}(x)\mu(x)\right)_1$','interpreter','latex','fontsize',16)
set(gca,'fontsize',16)
subplot(2,2,2);
surf(XX,YY,reshape(x_dot_at_x0(1,:),9,9))
xlabel('$x_1$','interpreter','latex','fontsize',16)
ylabel('$x_2$','interpreter','latex','fontsize',16)
zlabel('$\left(f(x)+ g(x)\mu(x)\right)_1$','interpreter','latex','fontsize',16)
set(gca,'fontsize',16)
subplot(2,2,3);
surf(XX,YY,reshape(x_dot_hat_at_x0(2,:),9,9))
xlabel('$x_1$','interpreter','latex','fontsize',16)
ylabel('$x_2$','interpreter','latex','fontsize',16)
zlabel('$\left(\hat{f}(x) + \hat{g}(x)\mu(x)\right)_2$','interpreter','latex','fontsize',16)
set(gca,'fontsize',16)
subplot(2,2,4);
surf(XX,YY,reshape(x_dot_at_x0(2,:),9,9))
xlabel('$x_1$','interpreter','latex','fontsize',16)
ylabel('$x_2$','interpreter','latex','fontsize',16)
zlabel('$\left(f(x) + g(x)\mu(x)\right)_2$','interpreter','latex','fontsize',16)
set(gca,'fontsize',16)

end

%% auxiliary functions
function out = oddLength(dt,tf)
    out = 0:dt:tf;
    if mod(numel(out),2)==0
        out = out(1:end-1);
    end
end
