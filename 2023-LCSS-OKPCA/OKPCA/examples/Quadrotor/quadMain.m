% This script uses occupation kernel principal component analysis for fault
% detection in an academic example.
%
% © Rushikesh Kamalapurkar
%
clear all; close all; clc;
addpath('../../lib')
%% Initialization
% Nominal system
n = 12; % State dimension
nI = 5; % Integral states for PID control
% Dataset parameters
tspan = 0:0.2:15.2; % Time span
M = 100; % # Training trajectories
MNormalTest = 20; % Number of normal test trajectories 
MFaultyTest = 20; % Number of faulty test trajectories 
noiseStandardDeviation = 0; % Standard deviation of added noise

% Kernel parameters
mu = 1000; % Kernel width
k = KernelRKHS('Gaussian',mu); % Kernel

% PCA parameters
N = 50; % Number of eigenvectors

% Fault detection parameters
thresholdMultiplier = 10; % Threshold = max training error times this

trials = 1; % Number of trials

% Matrices to store results
RTEST = zeros(MNormalTest+MFaultyTest,trials,numel(M));
trainingX0s = zeros(n+nI,M,trials);
testingNormalX0s = zeros(n+nI,MNormalTest,trials);
testingFaultyX0s = zeros(n+nI,MFaultyTest,trials);
RTRAIN = zeros(M,trials);

%% Repeated trials
for trial = 1:trials
    % Initial states for training data
    trainingX0 = [-1+2*rand(n,M);zeros(nI,M)];
    % Generate training data
    trainingPaths = zeros(n,length(tspan),length(trainingX0));
    for i = 1:M
        [~,temp] = ode45(@(t,z) normalModel(t,z),tspan,trainingX0(:,i));
        trainingPaths(:,:,i)=temp(:,1:n).';
    end
    tTraining = repmat(tspan.',1,M);
    
    % Initial states for normal test data
    testingNormalX0 = [-1+2*rand(n,MNormalTest);zeros(nI,MNormalTest)];
    % Generate normal test data
    testingNormalPaths = zeros(n,length(tspan),MNormalTest);
    for i = 1:MNormalTest
        [~,temp] = ode45(@(t,z) normalModel(t,z),tspan,testingNormalX0(:,i));
        testingNormalPaths(:,:,i)=temp(:,1:n).';
    end
    tNormalTest = repmat(tspan.',1,MNormalTest);
    
    % Initial states for abnormal test data
    testingFaultyX0 = [-1+2*rand(n,MFaultyTest);zeros(nI,MFaultyTest)];
    % Generate abnormal test data
    testingFaultyPaths = zeros(n,length(tspan),MFaultyTest);
    for i = 1:MFaultyTest
        [~,temp] = ode45(@(t,z) faultyModel(t,z),tspan,testingFaultyX0(:,i));
        testingFaultyPaths(:,:,i)=temp(:,1:n).';
    end
    tFaultyTest = repmat(tspan.',1,MFaultyTest);
    
    % Add noise
    testingNormalPaths = testingNormalPaths + noiseStandardDeviation*randn(size(testingNormalPaths));
    testingFaultyPaths = testingFaultyPaths + noiseStandardDeviation*randn(size(testingFaultyPaths));
    trainingPaths = trainingPaths + noiseStandardDeviation*randn(size(trainingPaths));

    % Store initial conditions
    trainingX0s(:,:,trial) = trainingX0;
    testingNormalX0s(:,:,trial) = testingNormalX0;
    testingFaultyX0s(:,:,trial) = testingFaultyX0;

    % OKPCA Reconstruction Error
    [RTest,RTrain] = OKPCAReconstructionError(k,trainingPaths,tTraining,cat(3,testingNormalPaths,testingFaultyPaths),cat(2,tNormalTest,tFaultyTest),N);

    % Compute and store 
    RTEST(:,trial) = RTest;
    RTRAIN(:,trial) = RTrain;
end

%% Fault detection
% Threshold for fault detection for each trial
epsilon = thresholdMultiplier*max(RTRAIN);

% False positives per trial
falsePositives = sum(RTEST(1:MNormalTest,:) > epsilon);

% False negatives per trial
falseNegatives = sum(RTEST(MNormalTest+1:MNormalTest+MFaultyTest,:) < epsilon);

%Total errors 
disp(['Total Errors = ' num2str(sum(falsePositives + falseNegatives))])

%% Plots
[~,bestTrial] = min(falsePositives + falseNegatives);
scatter(1:MNormalTest,log(RTEST(1:MNormalTest,bestTrial)),'b','filled');
hold on
scatter(1:MNormalTest,log(RTEST(MNormalTest+1:MNormalTest+MFaultyTest,bestTrial)),'r','filled');
line([1,MNormalTest],[log(epsilon(bestTrial)),log(epsilon(bestTrial))],'color','g','linewidth',2);
set(gca,'fontsize',14);
legend("Normal","Faulty","Threshold",'interpreter','latex','fontsize',14);
xlabel("Test point",'interpreter','latex','fontsize',14);
ylabel("Log reconstruction error",'interpreter','latex','fontsize',14);
xlim([1 20]);
hold off

figure
[~,worstTrial] = max(falsePositives + falseNegatives);
scatter(1:MNormalTest,log(RTEST(1:MNormalTest,worstTrial)),'b','filled');
hold on
scatter(1:MNormalTest,log(RTEST(MNormalTest+1:MNormalTest+MFaultyTest,worstTrial)),'r','filled');
line([1,MNormalTest],[log(epsilon(worstTrial)),log(epsilon(worstTrial))],'color','g','linewidth',2);
set(gca,'fontsize',14);
legend("Normal","Faulty","Threshold",'interpreter','latex','fontsize',14);
xlabel("Test point",'interpreter','latex','fontsize',14);
ylabel("Log reconstruction error",'interpreter','latex','fontsize',14);
xlim([1 20]);
hold off
