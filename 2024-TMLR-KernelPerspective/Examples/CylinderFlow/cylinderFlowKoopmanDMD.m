% This script uses finite element snapshots from the DMD book for fluid 
% flow past a cylinder and uses Koopman DMD to generate a predictive
% model from the data.
%
% © Rushikesh Kamalapurkar and Joel Rosenfeld
%
function cylinderFlowKoopmanDMD()
%% Set up paths
DATAPATH = '../../DATA';
addpath('../../lib')

%% Kernel selection
mu = 1;
K = KernelRKHS('Gaussian',mu); 

%% Load and format data for DMD
load([DATAPATH '/FLUIDS/CYLINDER_ALL.mat']);
DATA = VORTALL/max(vecnorm(VORTALL)); % Using normalized vorticity data
Length = 30; % Number of snapshots used for DMD. Should be between 1 and
             % 150, cannot be 150 if PlotPrediction is ON.
h = 0.02; % Time step
% Input values X
X = DATA(:,1:Length);
% Output values Y = F(X)
Y = DATA(:,2:Length+1);

%% Kernel DMD
[~,~,~,~,dr,~] = KoopmanDMD(X,Y,K,h,l);
[~,~,~,~,drW,~] = WilliamsKDMD(X,Y,K,h);

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
legend("Gonzalez et al.","Williams et al.");
xlabel("t [s]");
ylabel("Norm of the relative reconstruction error")
reconstructionError = [time reconstructionError];
save('KoopmanError.dat','reconstructionError','-ascii');
reconstructionErrorWilliams = [time reconstructionErrorWilliams];
save('WilliamsError.dat','reconstructionErrorWilliams','-ascii');
end