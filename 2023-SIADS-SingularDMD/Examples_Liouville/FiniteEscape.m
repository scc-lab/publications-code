close all
clear all

addpath('../Functions');

%% Generate data for DMD
xrange = 1;
numTraj = 10;
x0 = linspace(-xrange,xrange,numTraj).';
deltaT = 0.02;
tFinal = 0.5*ones(1,numTraj);
maxLength = length(0:deltaT:max(tFinal));
SampleTime = cell2mat(cellfun(@(x) [x;NaN(maxLength-length(x),1)],...
    arrayfun(@(x) (oddLength(deltaT,x)).',tFinal,'UniformOutput',false), 'UniformOutput', false));

xDot = @(t,x) 1+x.^2;
Trajectories=zeros(1,size(SampleTime,1),numTraj);

for i=1:numTraj
    T_i = SampleTime(~isnan(SampleTime(:,i)),i);
    [~,y]=ode45(@(t,x) xDot(t,x),T_i,x0(i));
    Trajectories(:,~isnan(SampleTime(:,i)),i)=y(:,1).';
end

GramMatrixRegularizationParameter = 1e-7;
%% Explicit kernel functions that can take 3D arrays as inputs (R2020b or
% newer version of MATLAB needed)
K = Kernel('Exponential',4.71);

%% Convergent Liouville DMD Eigenfunction Approach
[Modes,Eigenvalues,EigenfunctionCoeffs,ReconstructionFunction,VectorField] = ...
    LiouvilleDMD(K,Trajectories,SampleTime,1,GramMatrixRegularizationParameter);

%% Reconstruction
T = 0:0.01:0.7;
x = -0.5;
% Actual trajectory
[~,y]=ode45(@(t,x) xDot(t,x),T,x);
yr = zeros(size(y(:,1)));
for i=1:numel(T)
    yr(i,:) = real(ReconstructionFunction(T(i),x));
end
plot(T,y(:,1),T,yr,'LineWidth',2);
legend('$x(t)$','$\sum_{m=1}^{10} \hat \xi_{m} e^{\lambda_m t} \hat \varphi_m(x(0))$','interpreter','latex','location','northwest');    
xlabel('Time (s)','interpreter','latex')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [9 3]);
set(gcf, 'PaperPosition', [0 0 9 3]);
set(gca,'FontSize',16);
filename = ['finite-escape-' K.type '-liouville-reg-' num2str(GramMatrixRegularizationParameter) '-reconstruction.pdf'];
saveas(gcf,filename);
%% Vector field
xrange = 2;
pointsPerDim = 100;%10;
x = linspace(-xrange,xrange,pointsPerDim).';
xDotActual = 1 + x.^2;
xDotDMD = zeros(size(x));
for i=1:numel(x)
    xDotDMD(i) = VectorField(x(i));
end
figure
plot(x,abs(xDotActual - xDotDMD),'LineWidth',2);
legend('$\left\vert f(x) - \sum_{m=1}^{10} \lambda_m \hat \xi_{m} \hat \varphi_m(x)\right\vert$','Interpreter','latex');   
xlabel('$x$','interpreter','latex')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [9 3]);
set(gcf, 'PaperPosition', [0 0 9 3]);
set(gca,'FontSize',16);
filename = ['finite-escape-' K.type '-liouville-reg-' num2str(GramMatrixRegularizationParameter) '-vectorfield.pdf'];
saveas(gcf,filename);
%% functions
function out = oddLength(dt,tf)
    out = 0:dt:tf;
    if mod(numel(out),2)==0
        out = out(1:end-1);
    end
end