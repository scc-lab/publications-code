%% Paper title "Safe Model-Based Reinforcement Learning for Systems with Parametric Uncertainties" 
%% arXiv:2007.12666
%% Coded by S M Nahid Mahmud, MS Grad Student, Oklahoma State University.
%% nahid.mahmud@okstate.edu
%% Two-state dynamical system (Section 7.1 of the journal paper)

clear all
clc
ADP = load('Tranformed_States_Trajectory_RL_BF.mat');
GPOPS= load('Transformed_State_Trajectory_Jpops.mat');
FCL4 = load('States_RL_BF_FCL.mat');


% Phase portrait of the s-coordinates using BT MBRL and GPOPS II
figure(1)
p=plot(ADP.s(:,1),ADP.s(:,2),GPOPS.output.result.solution.phase.state(:,1),GPOPS.output.result.solution.phase.state(:,2),'Linewidth',1);
mrk1={'s','v','o','*','^','x'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('$s_{1}(t)$','interpreter','latex')
ylabel('$s_{2}(t)$','interpreter','latex')
xlim([-3.5 0])
set(gca,'FontSize',15)
l=legend('BT MBRL with optimal weights','GPOPS II','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 4]);
set(gcf, 'PaperPosition', [0 0 8 4]);
grid on 
% print -dpdf Time_Vs_TStates_RL_BF_FCL_Optimal.png
% print -dpdf Time_Vs_TStates_RL_BF_FCL_Optimal.pdf

% State trajectories of the c-coordinates using BT MBRL 
figure(2)

plot(FCL4.z(:,50),FCL4.z(:,51),'Linewidth',1)
rectangle('position',[-7 -5 12 12])
axis([-8 7 -7 8])
xlabel('$x_{1}(t)$','interpreter','latex')
ylabel('$x_{2}(t)$','interpreter','latex')
set(gca,'FontSize',15)
% l=legend('Original State Trajectory','interpreter','latex')
% set(l,'Interpreter','latex','Location','northeast','NumColumns',2);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 4]);
set(gcf, 'PaperPosition', [0 0 8 4]);
grid on 
% print -dpdf Main_States_RL_BF_FCL.png
% print -dpdf Main_States_RL_BF_FCL.pdf

% Trajectory of the ACtor-Critic weights 
figure(3);
p=plot(FCL4.t,FCL4.z(:,35),FCL4.t,FCL4.z(:,36),FCL4.t,FCL4.z(:,37),FCL4.t,FCL4.z(:,38),FCL4.t,FCL4.z(:,39),FCL4.t,FCL4.z(:,40),'Linewidth',1);
mrk1={'s','v','o','*','^','x'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('$t$','interpreter','latex')
ylabel({'$\hat{W}(t)$'},'interpreter','latex')
xlim([0 8])
ylim([0 10])
set(gca,'FontSize',15)
l=legend('$\hat{W_{a}}_{1}$','$\hat{W_{a}}_{2}$','$\hat{W_{a}}_{3}$','$\hat{W_{c}}_{1}$','$\hat{W_{c}}_{2}$','$\hat{W_{c}}_{3}$','interpreter','latex');
set(l,'Interpreter','latex','Location','northeast','NumColumns',2);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 4]);
set(gcf, 'PaperPosition', [0 0 8 4]);
grid on 
% print -dpdf Weights_Estimations_for_States_RL_BF_FCL.png
% print -dpdf Weights_Estimations_for_States_RL_BF_FCL.pdf

% Trajectory of the Parameter Estimation errors
figure(4);
p=plot(FCL4.t,(FCL4.P.theta(1)-FCL4.z(:,3)),FCL4.t,(FCL4.P.theta(2)-FCL4.z(:,4)),FCL4.t,(FCL4.P.theta(3)-FCL4.z(:,5)),FCL4.t,(FCL4.P.theta(4)-FCL4.z(:,6)),'Linewidth',1);
mrk1={'s','v','o','*','^','x'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('$t$','interpreter','latex')
ylabel({'$\tilde{\theta}(t)$'},'interpreter','latex')

set(gca,'FontSize',15)
l=legend('$\tilde{\theta}_{1}$','$\tilde{\theta}_{2}$','$\tilde{\theta}_{3}$','$\tilde{\theta}_{4}$','interpreter','latex');
set(l,'Interpreter','latex','Location','northeast','NumColumns',4);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 4]);
set(gcf, 'PaperPosition', [0 0 8 4]);
grid on 
% print -dpdf Parameter_Estimation_Errors_RL_BF_FCL.png
% print -dpdf Parameter_Estimation_Errors_RL_BF_FCL.pdf

% Trajectory of the Estimated parameters
figure(5);
p=plot(FCL4.t,FCL4.z(:,3),FCL4.t,FCL4.z(:,4),FCL4.t,FCL4.z(:,5),FCL4.t,FCL4.z(:,6),'Linewidth',1);
mrk1={'s','v','o','*','^','x'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('$t$','interpreter','latex')
ylabel({'$\hat{\theta}(t)$'},'interpreter','latex')
xlim ([0 4])
ylim ([-1 1])
set(gca,'FontSize',15)
hold on
p1=plot(FCL4.t,1*ones(size(FCL4.t)),'--','LineWidth',1);%creating a straight line at y = 1
set(p1,'Color',get(p(1),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(FCL4.t,-1*ones(size(FCL4.t)),'--','LineWidth',1);%creating a straight line at y = -1
set(p1,'Color',get(p(2),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(FCL4.t,-0.5*ones(size(FCL4.t)),'--','LineWidth',1);%creating a straight line at y = -0.5
set(p1,'Color',get(p(3),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(FCL4.t,0.5*ones(size(FCL4.t)),'--','LineWidth',1);%creating a straight line at y = 0.5
set(p1,'Color',get(p(4),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
l=legend('$\hat{\theta}_{1}$','$\hat{\theta}_{2}$','$\hat{\theta}_{3}$','$\hat{\theta}_{4}$','interpreter','latex');
set(l,'Interpreter','latex','Location','northeast','NumColumns',2);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 4]);
set(gcf, 'PaperPosition', [0 0 8 4]);
grid on 
% print -dpdf Parameter_Estimation_RL_BF_FCL.png
% print -dpdf Parameter_Estimation_RL_BF_FCL.pdf
