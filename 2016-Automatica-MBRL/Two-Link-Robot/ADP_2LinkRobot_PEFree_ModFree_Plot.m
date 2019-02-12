%% WcH

plot(WcH.time,WcH.signals.values,'LineWidth',2);
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{W}_{c}(t)$','FontSize',16,'Interpreter','latex')
title('Value Function Weights','FontSize',16)
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPEFree2LinkWc.pdf
close all

%% WaH

plot(WaH.time,WaH.signals.values,'LineWidth',2);
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{W}_{a}(t)$','FontSize',16,'Interpreter','latex')
title('Policy Weights','FontSize',16)
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPEFree2LinkWa.pdf
close all

%% Theta

p=plot(ThetaH.time,ThetaH.signals.values,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{\theta}(t)$','FontSize',14,'Interpreter','latex')
title('Unknown Parameters in Drift Dynamics','FontSize',16)
set(gca,'FontSize',16)
hold on
p1=plot(ThetaH.time,5.3*ones(size(ThetaH.time)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(ThetaH.time,1.1*ones(size(ThetaH.time)),'--','LineWidth',1);
set(p1,'Color',get(p(2),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(ThetaH.time,8.45*ones(size(ThetaH.time)),'--','LineWidth',1);
set(p1,'Color',get(p(3),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(ThetaH.time,2.35*ones(size(ThetaH.time)),'--','LineWidth',1);
set(p1,'Color',get(p(4),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
l=legend('$f_{d1}$','$f_{d2}$','$f_{s1}$','$f_{s2}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPEFree2LinkTheta.pdf
clear p
close all

%% State

p=plot(x.time,x.signals.values,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$x(t)$','FontSize',16,'Interpreter','latex')
title('State Trajectory','FontSize',16)
set(gca,'FontSize',16)
l=legend('$x_1$','$x_2$','$x_3$','$x_4$');
set(l,'Interpreter','latex','Location','southeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPEFree2LinkState.pdf
clear p
close all

%% Control

p=plot(u.time,u.signals.values,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$u(t)$','FontSize',16,'Interpreter','latex')
title('Control Trajectory','FontSize',16)
set(gca,'FontSize',16)
l=legend('$u_1$','$u_2$');
set(l,'Interpreter','latex','Location','southeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPEFree2LinkControl.pdf

%% Assumption 2

p=plot(cl1.time,cl1.signals.values,'LineWidth',2);
xlabel('Time (s)','FontSize',16)
ylabel('$\underline{y}(t)$','FontSize',16,'Interpreter','latex')
title('Satisfaction of Assumption 2','FontSize',16)
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPEFree2LinkSysIDSingularValue.pdf

%% Assumption 2

p=plot(cbar.time,cbar.signals.values,'LineWidth',2);
xlabel('Time (s)','FontSize',16)
ylabel('$\underline{c}(t)$','FontSize',16,'Interpreter','latex')
title('Satisfaction of Assumption 3','FontSize',16)
axis([0 100 0 10e-12])
set(gca,'FontSize',16)
p = gcf
set(p, 'PaperPositionMode', 'manual');
set(p, 'PaperUnits', 'inches');
set(p, 'PaperSize', [6 5]);
set(p, 'PaperPosition', [0 0 6 5]);
x_a = 0.2; y_a = 0.50; w_a = 0.3; h_a = 0.3;
ax = axes('Units', 'Normalized', ...
          'Position', [x_a, y_a, w_a, h_a], ...
          'Box', 'on', ...
          'LineWidth', 1, ...
          'Color', [1, 1, 1]);
hold on;
plot(cbar.time,cbar.signals.values);
xlabel('Time (s)','FontSize',10)
ylabel('$\underline{c}(t)$','FontSize',10,'Interpreter','latex')



set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPEFree2LinkBEExtraSingularValue.pdf