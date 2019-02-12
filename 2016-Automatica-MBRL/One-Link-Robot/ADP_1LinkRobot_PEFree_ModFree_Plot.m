Refsize = 3;
%% WcH
WcH=getElement(logsout,'WcH');
TT=WcH.Values.Time;
WW=WcH.Values.Data;
p=plot(TT,WW,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{W}_{c1}(t)$','FontSize',16,'Interpreter','latex')
title('Value Function Weights','FontSize',16)
set(gca,'FontSize',16)
hold on
p1=plot(TT,0.5*ones(size(TT)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',2)
hasbehavior(p1, 'legend', false);
p1=plot(TT,0*ones(size(TT)),'--','LineWidth',1);
set(p1,'Color',get(p(2),'Color'),'LineWidth',2)
hasbehavior(p1, 'legend', false);
p1=plot(TT,1*ones(size(TT)),'--','LineWidth',1);
set(p1,'Color',get(p(3),'Color'),'LineWidth',2)
hasbehavior(p1, 'legend', false);
l=legend('$\hat{W}_{c1,1}$','$\hat{W}_{c1,2}$','$\hat{W}_{c1,3}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPEFree1LinkWcH.pdf
clear p
close all

%% WaH
WaH=getElement(logsout,'WaH');
TT=WaH.Values.Time;
WW=WaH.Values.Data;
p=plot(TT,WW,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{W}_{c1}(t)$','FontSize',16,'Interpreter','latex')
title('Policy Weights','FontSize',16)
set(gca,'FontSize',16)
hold on
p1=plot(TT,0.5*ones(size(TT)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',2)
hasbehavior(p1, 'legend', false);
p1=plot(TT,0*ones(size(TT)),'--','LineWidth',1);
set(p1,'Color',get(p(2),'Color'),'LineWidth',2)
hasbehavior(p1, 'legend', false);
p1=plot(TT,1*ones(size(TT)),'--','LineWidth',1);
set(p1,'Color',get(p(3),'Color'),'LineWidth',2)
hasbehavior(p1, 'legend', false);
l=legend('$\hat{W}_{a1,1}$','$\hat{W}_{a1,2}$','$\hat{W}_{a1,3}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPEFree1LinkWaH.pdf
clear p
close all

%% Theta 
ThetaH=getElement(logsout,'SIDThetaHat');
TT=ThetaH.Values.Time;
WW=ThetaH.Values.Data;
p=plot(TT,WW,'LineWidth',2);
mrk1={'s','v','o','*','^','+'};
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
p1=plot(TT,-1*ones(size(TT)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',2)
hasbehavior(p1, 'legend', false);
p1=plot(TT,1*ones(size(TT)),'--','LineWidth',1);
set(p1,'Color',get(p(2),'Color'),'LineWidth',2)
hasbehavior(p1, 'legend', false);
p1=plot(TT,-0.5*ones(size(TT)),'--','LineWidth',1);
set(p1,'Color',get(p(3),'Color'),'LineWidth',2)
hasbehavior(p1, 'legend', false);
l=legend('$\hat{\theta}_1$','$\hat{\theta}_2$','$\hat{\theta}_3$','$\hat{\theta}_4$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPEFree1LinkTheta.pdf
clear p
close all

%% State
x=getElement(logsout,'x');
TT=x.Values.Time;
WW=x.Values.Data;
p=plot(TT,WW,'LineWidth',2);
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
l=legend('$x_1$','$x_2$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPEFree1LinkState.pdf
clear p
close all

%% Control
u=getElement(logsout,'u');
TT=u.Values.Time;
WW=u.Values.Data;
plot(TT,WW,'LineWidth',2)
xlabel('Time (s)','FontSize',16)
title('Control Trajectory','FontSize',16)
ylabel('$u\left(t\right)$','FontSize',16,'Interpreter','latex')
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPEFree1LinkControl.pdf
clear p
close all

%% SysIDSingular
cl1=getElement(logsout,'SIDRank');
TT=cl1.Values.Time;
WW=cl1.Values.Data;
plot(TT,WW,'LineWidth',2)
xlabel('Time (s)','FontSize',16)
title('Satisfaction of Assumption 2','FontSize',16)
ylabel('$\underline{y}(t)$','FontSize',16,'Interpreter','latex')
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPEFree1LinkSysIDSing.pdf
clear p
close all

%% ADPSingular
cbar=getElement(logsout,'cbar');
TT=cbar.Values.Time;
WW=cbar.Values.Data;
% p=figure;
% hold on
plot(TT,WW,'LineWidth',2)
% hold off
xlabel('Time (s)','FontSize',16)
title('Satisfaction of Assumption 3','FontSize',16)
ylabel('$\underline{c}(t)$','FontSize',16,'Interpreter','latex')
set(gca,'FontSize',16)
ylim([0 0.1]);
p = gcf
set(p, 'PaperPositionMode', 'manual');
set(p, 'PaperUnits', 'inches');
set(p, 'PaperSize', [6 5]);
set(p, 'PaperPosition', [0 0 6 5]);
x_a = 0.58; y_a = 0.58; w_a = 0.3; h_a = 0.3;
ax = axes('Units', 'Normalized', ...
          'Position', [x_a, y_a, w_a, h_a], ...
          'Box', 'on', ...
          'LineWidth', 1, ...
          'Color', [1, 1, 1]);
hold on;
plot(TT, WW);
xlabel('Time (s)','FontSize',10)
ylabel('$\underline{c}(t)$','FontSize',10,'Interpreter','latex')
print -dpdf ADPPEFree1LinkADPSing.pdf
clear p
close all


