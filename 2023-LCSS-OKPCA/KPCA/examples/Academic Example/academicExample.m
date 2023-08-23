% This script uses kernel principal component analysis for fault
% detection in an academic example. Hoffman's implementation of KPCA is
% used with minimal modification.
%
% © Rushikesh Kamalapurkar
%
clear all; close all; clc;
addpath('../../lib')

%% Initialization
dyn1 = @(t,x) -x + [x(2)*sin(pi/2*x(1)); x(1)*cos(pi/2*x(2))];
dyn2 = @(t,x) -x + [0.9*x(2)*sin(pi/5*x(1)); 0.8*x(1)*cos(pi/3*x(2))];
n = 2; % State dimension
tspan = 0:0.01:2; % Time span
M = 100; % # Training trajectories
MFaultyTest = 20; % # Faulty test trajectories 
MNormalTest = 20; % # Normal test trajectories
mu = 5; % Kernel width
N = 20; % Number of eigenvectors
noiseStandardDeviation = 0; % Noise standard deviation
trials = 1; % Number of repeated trials

% Storage matrices
falsePositives = zeros(trials,numel(M));
falseNegatives = zeros(trials,numel(M));
testingFaultyInitialParam = zeros(trials,MFaultyTest);
trainingInitialParams = zeros(trials,M);
RTest = zeros(trials,MNormalTest+MFaultyTest);
RTrain = zeros(trials,M);

%% Monte-carlo trials
for trial = 1:trials
    % Training data
    trainingInitialParam = 2*pi*rand(1,M);
    trainingX0 = [sin(trainingInitialParam);cos(trainingInitialParam)];
    % Generate training data
    trainingPaths = zeros(length(trainingX0),n*length(tspan));
    for i = 1:length(trainingX0)
        [~,temp] = ode45(dyn1,tspan,trainingX0(:,i));
        temp = temp.';
        trainingPaths(i,:)=temp(:)';
    end
    trainingPaths = trainingPaths + noiseStandardDeviation*randn(size(trainingPaths));
    
    % Faulty test data
    testingFaultyInitialParam = 2*pi*rand(1,MFaultyTest);
    testingFaultyX0 = [sin(testingFaultyInitialParam);cos(testingFaultyInitialParam)];
    % Generate faulty test data
    testingFaultyPaths = zeros(length(testingFaultyX0),n*length(tspan));
    for i = 1:length(testingFaultyX0)
        [~,temp] = ode45(dyn2,tspan,testingFaultyX0(:,i));
        temp = temp.';
        testingFaultyPaths(i,:)=temp(:)';
    end
    testingFaultyPaths = testingFaultyPaths + noiseStandardDeviation*randn(size(testingFaultyPaths));
    
    % Normal test data
    testingNormalInitialParam = 2*pi*rand(1,MNormalTest);
    testingNormalX0 = [sin(testingNormalInitialParam);cos(testingNormalInitialParam)];
    % Generate Normal test data
    testingNormalPaths = zeros(length(testingNormalX0),n*length(tspan));
    for i = 1:length(testingNormalX0)
        [~,temp] = ode45(dyn1,tspan,testingNormalX0(:,i));
        temp = temp.';
        testingNormalPaths(i,:)=temp(:)';
    end
    testingNormalPaths = testingNormalPaths + noiseStandardDeviation*randn(size(testingNormalPaths));

    trainingInitialParams(trial,:) = trainingInitialParam;
    testingFaultyInitialParam(trial,:) = testingFaultyInitialParam;
    testingNormalInitialParam(trial,:) = testingNormalInitialParam;
    
    % KPCA
    [rtrain,rtest] = kpcabound(trainingPaths,mu,N,[testingNormalPaths;testingFaultyPaths]);
    RTrain(trial,:) = rtrain;
    RTest(trial,:) = rtest;
    epsilon = 2*max(rtrain);
    falsePositives(trial,1) = nnz(rtest(1:MNormalTest) > epsilon);
    falseNegatives(trial,1) = nnz(rtest(MNormalTest+1:end) < epsilon);
end
scatter(1:MNormalTest,log(RTest(1,1:MNormalTest)),'b','filled')
hold on
scatter(1:MNormalTest,log(RTest(1,MNormalTest+1:MNormalTest+MFaultyTest)),'r','filled')
line([1,20],log([2*max(RTrain(1,:)),2*max(RTrain(1,:))]),'color','g','linewidth',2)
