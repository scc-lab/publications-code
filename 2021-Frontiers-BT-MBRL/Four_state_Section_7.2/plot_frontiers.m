%% Paper title "Safe Model-Based Reinforcement Learning for Systems with Parametric Uncertainties" 
%% arXiv:2007.12666
%% Coded by S M Nahid Mahmud, MS Grad Student, Oklahoma State University.
%% nahid.mahmud@okstate.edu
%% Four-state dynamical system (Section 7.2 of the journal paper)

clear all
clc
ADP = load('States_RL_BF_manipulator.mat');
GPOPS= load('Time_vs_States_Jpops.mat');
FCL4 = load('States_RL_BF_FCL_manipulator_final_last.mat');



% State trajectories of the s-coordinates using BT MBRL and GPOPS II
figure(1)
p=plot(ADP.t,ADP.z(:,1:2),GPOPS.time,GPOPS.state(:,1:2),'Linewidth',1);
mrk1={'s','v','o','*','^','x','+','.','d'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('$t$','interpreter','latex')
ylabel({'$s(t)$'},'interpreter','latex')
xlim([0 150])
ylim([-2 1])
set(gca,'FontSize',15)
l=legend('BT MBRL, $s_{1}$','BT MBRL, $s_{2}$','GPOPS II, $s_{1}$',...
    'GPOPS II, $s_{2}$','interpreter','latex')
set(l,'Interpreter','latex','Location','southeast','NumColumns',2);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 2]);
set(gcf, 'PaperPosition', [0 0 8 2]);
grid on 
% print -dpdf Time_Vs_TStates_RL_BF_FCL_Optimal_robot1.png
 print -dpdf Time_Vs_TStates_RL_BF_FCL_Optimal_robot1.pdf

 
figure(2)
p=plot(ADP.t,ADP.z(:,3:4),GPOPS.time,GPOPS.state(:,3:4),'Linewidth',1);
mrk1={'s','v','o','*','^','x','+','.','d'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,50,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('$t$','interpreter','latex')
ylabel({'$s(t)$'},'interpreter','latex')
xlim([0 20])
ylim([-1 2])
set(gca,'FontSize',15)
l=legend('BT MBRL, $s_{3}$','BT MBRL, $s_{4}$',...
    'GPOPS II, $s_{3}$','GPOPS II, $s_{4}$','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast','NumColumns',2);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 2]);
set(gcf, 'PaperPosition', [0 0 8 2]);
grid on 
% print -dpdf Time_Vs_TStates_RL_BF_FCL_Optimal_robot2.png
%  print -dpdf Time_Vs_TStates_RL_BF_FCL_Optimal_robot2.pdf
 
 
 
 
 
 
 
% State trajectories of the x-coordinates
figure(3);
p=plot(FCL4.t,FCL4.z(:,FCL4.P.n+FCL4.P.l+FCL4.P.n*FCL4.P.l+FCL4.P.l^2+FCL4.P.l+2*FCL4.P.W+FCL4.P.G^2+1:FCL4.P.n+FCL4.P.l+FCL4.P.n*FCL4.P.l+FCL4.P.l^2+FCL4.P.l+2*FCL4.P.W+FCL4.P.G^2+FCL4.P.n),'Linewidth',1);
mrk1={'s','v','o','*','^','x'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,20,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('$t$','interpreter','latex')
ylabel({'$x(t)$'},'interpreter','latex')
xlim([0 65])
ylim([-8 8])
set(gca,'FontSize',15)
hold on
p1=plot(FCL4.t,-7*ones(size(FCL4.t)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(FCL4.t,5*ones(size(FCL4.t)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(FCL4.t,-5*ones(size(FCL4.t)),'--','LineWidth',1);
set(p1,'Color',get(p(3),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(FCL4.t,7*ones(size(FCL4.t)),'--','LineWidth',1);
set(p1,'Color',get(p(3),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
l=legend('$x_{1}$','$x_{2}$','$x_{3}$','$x_{4}$','interpreter','latex');
set(l,'Interpreter','latex','Location','northeast','NumColumns',2);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 4]);
set(gcf, 'PaperPosition', [0 0 8 4]);
grid on 
% 
% print -dpdf State_RL_BF_FCL_robot.png
%  print -dpdf State_RL_BF_FCL_robot.pdf
% 
% 

% Trajectory of the Critic weights 
figure(4);
p=plot(FCL4.t,FCL4.z(:,FCL4.P.n+FCL4.P.l+FCL4.P.n*FCL4.P.l+FCL4.P.l^2+FCL4.P.l+FCL4.P.W+1:FCL4.P.n+FCL4.P.l+FCL4.P.n*FCL4.P.l+FCL4.P.l^2+FCL4.P.l+2*FCL4.P.W),'Linewidth',1);
mrk1={'s','v','o','*','^','x','+','.','p','d'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,20,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('$t$','interpreter','latex')
ylabel({'$\hat{W_{c}}(t)$'},'interpreter','latex')
xlim([0 12])
set(gca,'FontSize',15)
l=legend('$\hat{W_{c}}_{1}$','$\hat{W_{c}}_{2}$','$\hat{W_{c}}_{3}$','$\hat{W_{c}}_{4}$','$\hat{W_{c}}_{5}$',...
    '$\hat{W_{c}}_{6}$','$\hat{W_{c}}_{7}$','$\hat{W_{c}}_{8}$','$\hat{W_{c}}_{9}$','$\hat{W_{c}}_{10}$','interpreter','latex');
set(l,'Interpreter','latex','Location','Best','NumColumns',5);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 4]);
set(gcf, 'PaperPosition', [0 0 8 4]);
grid on 
% print -dpdf Weights_Estimations_for_States_RL_BF_FCL_robot.png
% print -dpdf Weights_Estimations_for_States_RL_BF_FCL_robot.pdf
% % 
% % figure(4);
% % p=plot(FCL4.t,FCL4.z(:,FCL4.P.n+FCL4.P.l+FCL4.P.n*FCL4.P.l+FCL4.P.l^2+FCL4.P.l+1:FCL4.P.n+FCL4.P.l+FCL4.P.n*FCL4.P.l+FCL4.P.l^2+FCL4.P.l+FCL4.P.W)-...
% %     FCL4.z(:,FCL4.P.n+FCL4.P.l+FCL4.P.n*FCL4.P.l+FCL4.P.l^2+FCL4.P.l+FCL4.P.W+1:FCL4.P.n+FCL4.P.l+FCL4.P.n*FCL4.P.l+FCL4.P.l^2+FCL4.P.l+2*FCL4.P.W),'Linewidth',1);
% % mrk1={'s','v','o','*','^','x','+','.','p','d'};
% % mrk=(mrk1(1,1:size(p,1)))';
% % set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5;
% % f_nummarkers(p,20,1);
% % for i=1:size(p,1)
% %     hasbehavior(p(i), 'legend', false);
% % end
% % xlabel('Time, t(s)','interpreter','latex')
% % ylabel({'$\tilde{W}(t)$'},'interpreter','latex')
% % set(gca,'FontSize',10)
% % l=legend('$\tilde{W}_{1}$','$\tilde{W}_{2}$','$\tilde{W}_{3}$','$\tilde{W}_{4}$','$\tilde{W}_{5}$','$\tilde{W}_{6}$',...
% %     '$\tilde{W}_{7}$','$\tilde{W}_{8}$','$\tilde{W}_{9}$','$\tilde{W}_{10}$','interpreter','latex');
% % set(l,'Interpreter','latex','Location','northeast','NumColumns',5);
% % set(gcf, 'PaperPositionMode', 'manual');
% % set(gcf, 'PaperUnits', 'inches');
% % set(gcf, 'PaperSize', [4 2]);
% % set(gcf, 'PaperPosition', [0 0 4 2]);
% % % grid on 
% % print -dpdf Weights_Estimations_Errors_RL_BF_FCL_robot.png
% % print -dpdf Weights_Estimations_Errors_RL_BF_FCL_robot.pdf
% 


% Trajectory of the Estimated parameters
figure(5);
p=plot(FCL4.t,FCL4.z(:,FCL4.P.n+1:FCL4.P.n+FCL4.P.l),'Linewidth',1);
mrk1={'s','v','o','*','^','x'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,20,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('$t$','interpreter','latex')
ylabel({'$\hat{\theta}(t)$'},'interpreter','latex')
xlim([0 16 ])
set(gca,'FontSize',15)
hold on
p1=plot(FCL4.t,5.3*ones(size(FCL4.t)),'--','LineWidth',1); %creating a straight line at y = 5.3
set(p1,'Color',get(p(1),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(FCL4.t,1.1*ones(size(FCL4.t)),'--','LineWidth',1); %creating a straight line at y = 1.1
set(p1,'Color',get(p(2),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(FCL4.t,8.45*ones(size(FCL4.t)),'--','LineWidth',1); %creating a straight line at y = 8.45
set(p1,'Color',get(p(3),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(FCL4.t,2.35*ones(size(FCL4.t)),'--','LineWidth',1); %creating a straight line at y = 2.35
set(p1,'Color',get(p(4),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
l=legend('$\hat{\theta}_{1}$','$\hat{\theta}_{2}$','$\hat{\theta}_{3}$','$\hat{\theta}_{4}$','interpreter','latex');
set(l,'Interpreter','latex','Location','northeast','NumColumns',2);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 4]);
set(gcf, 'PaperPosition', [0 0 8 4]);
grid on 
% % print -dpdf Parameter_Estimation_RL_BF_FCL_robot.png
%  print -dpdf Parameter_Estimation_RL_BF_FCL_robot.pdf
