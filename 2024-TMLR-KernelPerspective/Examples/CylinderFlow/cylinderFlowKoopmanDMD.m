% This script uses finite element snapshots from the DMD book for fluid 
% flow past a cylinder and uses Koopman DMD to generate a predictive
% model from the data.
%
% © Rushikesh Kamalapurkar and Joel Rosenfeld
%
function cylinderFlowKoopmanDMD()
	
DATAPATH = '../../DATA';
addpath('../../lib')

%% Kernel selection
% mu = 5;
% l = 25e-12; % Gram matrix regularization parameter
mu = 1;
l = 0; % Gram matrix regularization parameter
K = KernelRKHS('Gaussian',mu); 

%% Load and format data for DMD
load([DATAPATH '/FLUIDS/CYLINDER_ALL.mat']);
DATA = VORTALL/max(vecnorm(VORTALL)); % Using normalized vorticity data to compare with the DMD book
Length = 30; % Number of snapshots used for DMD. Should be between 1 and
             % 150, cannot be 150 if PlotPrediction is ON.
h = 0.02; % Time step
Width = 449; % Width needed for plotting only

Dimension = size(DATA,1);
% Input values X
X = DATA(:,1:Length);
% Output values Y = F(X)
Y = DATA(:,2:Length+1);

%% Kernel DMD
tic
[~,~,~,~,dr,~] = KoopmanDMD(X,Y,K,h,l);
[~,~,~,~,drW,~] = WilliamsKDMD(X,Y,K,h);
toc
%% Plots for the paper
reconstructionError = zeros(size(DATA,2)-1,1);
reconstructionErrorWilliams = zeros(size(DATA,2)-1,1);
time = zeros(size(DATA,2)-1,1);
x = DATA(:,1);
xW = x;
for i=1:numel(reconstructionError)
    time(i) = h*i;
    reconstructionError(i) = norm(DATA(:,i+1) - dr(i,x));
    reconstructionErrorWilliams(i) = norm(DATA(:,i+1) - drW(i,xW));
end
plot(reconstructionError);hold on;plot(reconstructionErrorWilliams);hold off;
reconstructionError = [time reconstructionError];
save('KoopmanError.dat','reconstructionError','-ascii');
reconstructionErrorWilliams = [time reconstructionErrorWilliams];
save('WilliamsError.dat','reconstructionErrorWilliams','-ascii');
end