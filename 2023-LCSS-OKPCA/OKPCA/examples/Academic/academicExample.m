% This script uses occupation kernel principal component analysis for fault
% detection in an academic example.
%
% © Rushikesh Kamalapurkar
%
clear all; close all; clc;
addpath('../../lib')
%% Initialization
% Nominal system
n = 2; % State dimension
dyn1 = @(t,x) -x + [x(2)*sin(pi/2*x(1)); x(1)*cos(pi/2*x(2))];

% Faulty system
dyn2 = @(t,x) -x + [0.9*x(2)*sin(pi/5*x(1)); 0.8*x(1)*cos(pi/3*x(2))];

% Dataset parameters
tspan = 0:0.01:2; % Time span
M = 100; % # Training trajectories
% M = [50,100,150,300]; % # Training trajectories
MNormalTest = 20; % Number of normal test trajectories 
MFaultyTest = 20; % Number of faulty test trajectories 
noiseStandardDeviation = 0; % Standard deviation of added noise

% Kernel parameters
mu = 0.6; % Kernel width
k = KernelRKHS('Gaussian',mu); % Kernel

% PCA parameters
N = 20; % Number of eigenvectors

% Fault detection parameters
thresholdMultiplier = 2; % Threshold = max training error times this

trials = 1; % Number of trials

% Matrices to store results
RTEST = zeros(MNormalTest+MFaultyTest,trials,numel(M));
testingNormalInitialParams = zeros(trials,MNormalTest,numel(M));
testingFaultyInitialParam = zeros(trials,MFaultyTest,numel(M));
if numel(M)==1
    RTRAIN = zeros(M,trials);
    trainingInitialParams = zeros(trials,M);
else
    RTRAIN = cell(numel(M),1);
    trainingInitialParams = cell(numel(M),1);
    for ii=1:numel(M)
        RTRAIN{ii,1} = zeros(M(ii),trials);
        trainingInitialParams{ii,1} = zeros(trials,M(ii));
    end
end

%% Repeated trials
for ii = 1:numel(M)
    for trial = 1:trials
        % Initial states for training data
        trainingInitialParam = 2*pi*rand(1,M(ii));
        trainingX0 = [sin(trainingInitialParam);cos(trainingInitialParam)];
        % Generate training data
        trainingPaths = zeros(n,length(tspan),length(trainingX0));
        for i = 1:length(trainingX0)
            [~,temp] = ode45(dyn1,tspan,trainingX0(:,i));
            trainingPaths(:,:,i)=temp';
        end
        tTraining = repmat(tspan.',1,M(ii));
        
        % Initial states for normal test data
        testingNormalInitialParam = 2*pi*rand(1,MNormalTest);
        testingNormalX0 = [sin(testingNormalInitialParam);cos(testingNormalInitialParam)];
        % Generate normal test data
        testingNormalPaths = zeros(n,length(tspan),length(testingNormalX0));
        for i = 1:length(testingNormalX0)
            [~,temp] = ode45(dyn1,tspan,testingNormalX0(:,i));
            testingNormalPaths(:,:,i)=temp';
        end
        tNormalTest = repmat(tspan.',1,MNormalTest);
        
        % Initial states for abnormal test data
        testingFaultyInitialParam = 2*pi*rand(1,MFaultyTest);
        testingFaultyX0 = [sin(testingFaultyInitialParam);cos(testingFaultyInitialParam)];
        % Generate abnormal test data
        testingFaultyPaths = zeros(n,length(tspan),length(testingFaultyX0));
        for i = 1:length(testingFaultyX0)
            [~,temp] = ode45(dyn2,tspan,testingFaultyX0(:,i));
            testingFaultyPaths(:,:,i)=temp';
        end
        tFaultyTest = repmat(tspan.',1,MFaultyTest);
        
        % Add noise
        trainingPaths = trainingPaths + noiseStandardDeviation*randn(size(trainingPaths));
        testingNormalPaths = testingNormalPaths + noiseStandardDeviation*randn(size(testingNormalPaths));
        testingFaultyPaths = testingFaultyPaths + noiseStandardDeviation*randn(size(testingFaultyPaths));
        
        % Store initial conditions
        if numel(M) == 1
            trainingInitialParams(trial,:) = trainingInitialParam;
        else
            trainingInitialParams{ii}(trial,:) = trainingInitialParam;
        end
        testingNormalInitialParam(trial,:,ii) = testingNormalInitialParam;
        testingFaultyInitialParam(trial,:,ii) = testingFaultyInitialParam;

        % OKPCA Reconstruction Error
        [RTest,RTrain] = OKPCAReconstructionError(k,trainingPaths,tTraining,cat(3,testingNormalPaths,testingFaultyPaths),cat(2,tNormalTest,tFaultyTest),N);
        
        % Store reconstruction errors
        RTEST(:,trial,ii) = RTest;
        if numel(M) == 1
            RTRAIN(:,trial) = RTrain;
        else
            RTRAIN{ii}(:,trial) = RTrain;
        end
    end
end

%% Fault detection
if numel(M) == 1
    % Threshold for fault detection for each trial
    epsilon = thresholdMultiplier*max(RTRAIN);
else
    epsilon = zeros(1,trials,numel(M));
    for ii=1:numel(M)
        epsilon(1,:,ii) = thresholdMultiplier*max(RTRAIN{ii});
    end
end

% False positives per trial
falsePositives = sum(RTEST(1:MNormalTest,:,:) > epsilon);

% False negatives per trial
falseNegatives = sum(RTEST(MNormalTest+1:MNormalTest+MFaultyTest,:,:) < epsilon);

%Total errors 
disp(['Total Errors = ' num2str(sum(falsePositives + falseNegatives,'all'))])

%% Plots
if numel(M) == 1
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
    
    figure
    hold on
    handle1 = plot(squeeze(trainingPaths(1,:,:)),squeeze(trainingPaths(2,:,:)),'g');
    handle2 = plot(squeeze(testingNormalPaths(1,:,:)),squeeze(testingNormalPaths(2,:,:)),'b');
    handle3 = plot(squeeze(testingFaultyPaths(1,:,:)),squeeze(testingFaultyPaths(2,:,:)),'r');
    set(gca,'fontsize',14);
    legend([handle1(1),handle2(1),handle3(1)],'Training Data','Normal Test Data','Faulty Test Data','interpreter','latex','fontsize',14)
    xlabel("$x_1$",'interpreter','latex','fontsize',14);
    ylabel("$x_2$",'interpreter','latex','fontsize',14);
    hold off
end