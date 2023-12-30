% This script generates trajectories of a two link robot and
% uses control Liouville DMD to generate a predictive model for a given 
% feedback controller.
%
% © Rushikesh Kamalapurkar and Zachary Morrison
%
function controlledTwoLinkRobotManipulatorKoopman()
rng(1) % for reproducibility
addpath('../../lib')
%% Generate Trajectories
n = 4; % Number of dimensions that f maps from/to
m = 2; % Dimensions of the controller
f = @(x) ...
    [x(3);
     x(4);
    (1/(0.196^2 - 3.473*0.196 + 0.242^2*cos(x(2))^2)) *...
    [-0.196                   0.196 + 0.242*cos(x(2));
      0.196+0.242*cos(x(2))   -3.473-2*0.242*cos(x(2))] *...
    ((-[-0.242*sin(x(2))*x(4)   -0.242*sin(x(2))*(x(3)+x(4));
         0.242*sin(x(2))*x(3)    0                           ] - ...
    diag([5.3,1.1])) *...
    [x(3);
     x(4)] - ...
    [8.45*tanh(x(3));
     2.35*tanh(x(4))])];
g = @(x) ...
    [0 0;
     0 0;
    (1/(0.196^2 - 3.473*0.196 + 0.242^2*cos(x(2))^2))*...
    [-0.196                     0.196 + 0.242*cos(x(2))  ;...
      0.196 + 0.242*cos(x(2))  -3.473 - 2*0.242*cos(x(2))]];
IV_selection = 'random'; 
samp_min = -1;
samp_max = 1;
M = 300;
if strcmp(IV_selection,'random')
    % Get TotalTrajectories random IV's.
    X = samp_min + (samp_max - samp_min)*rand(n, M);
elseif strcmp(IV_selection,'halton')
    % Get TotalTrajectories halton sequence
    haltonseq = @(n,d) net(haltonset(d),n);
    halton = haltonseq(M, n);
    X = samp_min + (samp_max - samp_min)*halton.';
else
    error('Unknown IV selection mode %s', IV_selection)
end
ts = 0.1;
U = -2+4*rand(m,M);
Y=zeros(size(X));
for i = 1:M
    F = @(x,u) f(x) + g(x) * u; % The update function
    [~,y] = ode45(@(t,x) F(x,U(:,i)),[0,ts],X(:,i));
    Y(:,i) = y(end,:).';
end

%% Kernels
kT = 10000;
k = 10000;
l = 1e-6;

K=KernelvvRKHS('Exponential',k*ones(m+1,1));
KT=KernelRKHS('Exponential',kT);
%% Feedback controller
mu = @(x) cat(1, -5*x(1,:,:) - 5*x(2,:,:), -15*x(1,:,:) - 15*x(2,:,:));

%% DCLDMD
[~,~,~,~,dr,fHat] = ControlKoopmanDMD(KT,K,X,U,Y,ts,mu,l);

%% Indirect discrete reconstruction
t_pred = 0:0.05:6;
y = zeros(n,numel(t_pred));
y(:,1) = [1;-1;1;-1];
y_pred = zeros(n,numel(t_pred));
y_pred(:,1) = [1;-1;1;-1];
for i=1:numel(t_pred)-1
y_pred(:,i+1) = dr(1,y_pred(:,i));
[~,temp] = ode45(@(t,x) f(x) + g(x) * mu(x),[0,ts],y(:,i));
y(:,i+1) = temp(end,:).';
end

% Plots
plot(t_pred,y,'linewidth',2)
hold on
set(gca,'ColorOrderIndex',1)
plot(t_pred,y_pred,'--','linewidth',2)
hold off
xlabel('Time (s)')
set(gca,'fontsize',16)
legend('$x_1(t)$','$x_2(t)$','$x_3(t)$','$x_4(t)$','$\hat{x}_1(t)$','$\hat{x}_2(t)$','$\hat{x}_3(t)$','$\hat{x}_4(t)$',...
'interpreter','latex','fontsize',16,'location','southeast')
% TLRobotDiscreteReconCompare = [t_pred',y_pred',y']; % data for tikzplot
% TLRobotDiscreteReconError = [t_pred',y_pred'-y']; 
% save('TLRobotDiscreteReconCompare.dat','TLRobotDiscreteReconCompare','-ascii');
% save('TLRobotDiscreteReconError.dat','TLRobotDiscreteReconError','-ascii');
end

%% auxiliary functions
function out = oddLength(dt,tf)
    out = 0:dt:tf;
    if mod(numel(out),2)==0
        out = out(1:end-1);
    end
end
