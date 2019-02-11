Refsize = 3;
%% WcH1
p=plot(WcH1.time,WcH1.signals.values,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{W}_{c1}(t)$','FontSize',16,'Interpreter','latex')
title('Player 1 Value Function Weights','FontSize',16)
set(gca,'FontSize',16)
hold on
p1=plot(WcH1.time,0.5*ones(size(WcH1.time)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
p1=plot(WcH1.time,0*ones(size(WcH1.time)),'--','LineWidth',1);
set(p1,'Color',get(p(2),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
p1=plot(WcH1.time,1*ones(size(WcH1.time)),'--','LineWidth',1);
set(p1,'Color',get(p(3),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
l=legend('$\hat{W}_{c1,1}$','$\hat{W}_{c1,2}$','$\hat{W}_{c1,3}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
if model==0
    print -dpdf ADPPEFree1LinkWc1Nash.pdf
else
    print -dpdf ADPPEFree1LinkWc1NashACC.pdf
end
clear p
close all

%% WcH2
p=plot(WcH2.time,WcH2.signals.values,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{W}_{c2}(t)$','FontSize',16,'Interpreter','latex')
title('Player 2 Value Function Weights','FontSize',16)
set(gca,'FontSize',16)
hold on
p1=plot(WcH2.time,0.25*ones(size(WcH2.time)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
p1=plot(WcH2.time,0*ones(size(WcH2.time)),'--','LineWidth',1);
set(p1,'Color',get(p(2),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
p1=plot(WcH2.time,0.5*ones(size(WcH2.time)),'--','LineWidth',1);
set(p1,'Color',get(p(3),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
l=legend('$\hat{W}_{c2,1}$','$\hat{W}_{c2,2}$','$\hat{W}_{c2,3}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
if model==0
    print -dpdf ADPPEFree1LinkWc2Nash.pdf
else
print -dpdf ADPPEFree1LinkWc2NashACC.pdf
end
clear p
close all

%% WaH1

p=plot(WaH1.time,WaH1.signals.values,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{W}_{a1}(t)$','FontSize',16,'Interpreter','latex')
title('Player 1 Actor Weights','FontSize',16)
set(gca,'FontSize',16)
hold on
p1=plot(WaH1.time,0.5*ones(size(WaH1.time)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
p1=plot(WaH1.time,0*ones(size(WaH1.time)),'--','LineWidth',1);
set(p1,'Color',get(p(2),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
p1=plot(WaH1.time,1*ones(size(WaH1.time)),'--','LineWidth',1);
set(p1,'Color',get(p(3),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
l=legend('$\hat{W}_{a1,1}$','$\hat{W}_{a1,2}$','$\hat{W}_{a1,3}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
if model==0
    print -dpdf ADPPEFree1LinkWa1Nash.pdf
else
print -dpdf ADPPEFree1LinkWa1NashACC.pdf
end
clear p
close all

%% WaH2

p=plot(WaH2.time,WaH2.signals.values,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{W}_{a2}(t)$','FontSize',16,'Interpreter','latex')
title('Player 2 Actor Weights','FontSize',16)
set(gca,'FontSize',16)
hold on
p1=plot(WaH2.time,0.25*ones(size(WaH2.time)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
p1=plot(WaH2.time,0*ones(size(WaH2.time)),'--','LineWidth',1);
set(p1,'Color',get(p(2),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
p1=plot(WaH2.time,0.5*ones(size(WaH2.time)),'--','LineWidth',1);
set(p1,'Color',get(p(3),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
l=legend('$\hat{W}_{a2,1}$','$\hat{W}_{a2,2}$','$\hat{W}_{a2,3}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
if model==0
    print -dpdf ADPPEFree1LinkWa2Nash.pdf
else
print -dpdf ADPPEFree1LinkWa2NashACC.pdf
end
clear p
close all

%% Theta
if model==0
p=plot(ThetaH.time,ThetaH.signals.values,'LineWidth',2);
mrk1={'s','v','o','*','^','d'};
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
p1=plot(ThetaH.time,1*ones(size(ThetaH.time)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
p1=plot(ThetaH.time,-2*ones(size(ThetaH.time)),'--','LineWidth',1);
set(p1,'Color',get(p(2),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
p1=plot(ThetaH.time,-0.5*ones(size(ThetaH.time)),'--','LineWidth',1);
set(p1,'Color',get(p(3),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
p1=plot(ThetaH.time,-1*ones(size(ThetaH.time)),'--','LineWidth',1);
set(p1,'Color',get(p(4),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
p1=plot(ThetaH.time,0.25*ones(size(ThetaH.time)),'--','LineWidth',1);
set(p1,'Color',get(p(5),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
p1=plot(ThetaH.time,0.25*ones(size(ThetaH.time)),'-.','LineWidth',1);
set(p1,'Color',get(p(6),'Color'),'LineWidth',Refsize)
hasbehavior(p1, 'legend', false);
l=legend('$\hat{\theta}_1$','$\hat{\theta}_2$','$\hat{\theta}_3$','$\hat{\theta}_4$','$\hat{\theta}_5$','$\hat{\theta}_6$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPEFree1LinkThetaNashACC.pdf
clear p
close all
end
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
l=legend('$x_1$','$x_2$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
if model==0
    print -dpdf ADPPEFree1LinkStateNash.pdf
else
print -dpdf ADPPEFree1LinkStateNashACC.pdf
end
clear p
close all

%% Control

p=plot(u1.time,u1.signals.values,u2.time,u2.signals.values,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{u}\left(t\right)$','FontSize',16,'Interpreter','latex')
title('Control Policies','FontSize',16)
l=legend('Player 1','Player 2');
set(l,'Location','southeast');
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
if model==0
    print -dpdf ADPPEFree1LinkControlNash.pdf
else
print -dpdf ADPPEFree1LinkControlNashACC.pdf
end
clear p
close all

%% State Error
p=plot(x.time,x.signals.values-xStar.signals.values,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$x(t)-x^*(t)$','FontSize',16,'Interpreter','latex')
title('State Trajectory Error','FontSize',16)
set(gca,'FontSize',16)
l=legend('$x_1-x^*_1$','$x_2-x^*_2$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
if model==0
    print -dpdf ADPPEFree1LinkStateErrorNash.pdf
else
print -dpdf ADPPEFree1LinkStateErrorNashACC.pdf
end
clear p
close all

%% Control Error
p=plot(u1.time,u1.signals.values-u1Star.signals.values,u2.time,u2.signals.values-u2Star.signals.values,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{u}\left(t\right)-u^*\left(t\right)$','FontSize',16,'Interpreter','latex')
title('Control Estimation Error','FontSize',16)
l=legend('Player 1','Player 2');
set(l,'Location','southeast');
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
if model==0
    print -dpdf ADPPEFree1LinkControlErrorNash.pdf
else
print -dpdf ADPPEFree1LinkControlErrorNashACC.pdf
end
clear p
close all
