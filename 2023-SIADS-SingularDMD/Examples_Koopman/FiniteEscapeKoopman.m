close all
clear all

addpath('../Functions');

%% Finite escape trajectory generator
xrange = 1;
numTraj = 15;
x0 = linspace(-xrange,xrange,numTraj).';
deltaT = 0.25;
tFinal = 0.5;
X=[];
Y=[];
T = 0:deltaT:tFinal;
for i=1:numTraj
    xDot = @(t,x) 1+x.^2;
    [t,y] = ode45(@(t,x) xDot(t,x),T,x0(i));
    X = [X y(1:end-1,:).'];
    Y = [Y y(2:end,:).'];
end

%% Liouville DMD
mu = 1.7; % 1.7 for Williams, 2.5 for Koopman
l=1e-5;
K = Kernel('Exponential',mu);
[~,~,~,ContinuousReconstruction,DiscreteReconstruction] = ...
    WilliamsKDMD(X,Y,K,deltaT,l);

%% Reconstruction
T = 0:0.01:0.7;
x = -0.5;
% Actual trajectory
[~,y]=ode45(@(t,x) xDot(t,x),T,x);
yr = zeros(size(y(:,1)));
for i=1:numel(T)
    yr(i,:) = real(ContinuousReconstruction(T(i),x));
end
plot(T,y(:,1),T,yr,'LineWidth',2);
legend('$x(t)$','$\hat{x}(t)$','interpreter','latex','location','northwest');    
xlabel('Time (s)','interpreter','latex')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [9 3]);
set(gcf, 'PaperPosition', [0 0 9 3]);
set(gca,'FontSize',16);
filename = ['finite-escape-' K.type '-williams-KDMD-reconstruction.pdf'];
saveas(gcf,filename);