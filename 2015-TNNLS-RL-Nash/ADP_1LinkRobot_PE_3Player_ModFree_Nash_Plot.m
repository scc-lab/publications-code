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
l=legend('$\hat{W}_{c1,1}$','$\hat{W}_{c1,2}$','$\hat{W}_{c1,3}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPE1LinkWc1Nash3Player.pdf
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
l=legend('$\hat{W}_{c2,1}$','$\hat{W}_{c2,2}$','$\hat{W}_{c2,3}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPE1LinkWc2Nash3Player.pdf
clear p
close all

%% WcH3
p=plot(WcH3.time,WcH3.signals.values,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{W}_{c3}(t)$','FontSize',16,'Interpreter','latex')
title('Player 3 Value Function Weights','FontSize',16)
set(gca,'FontSize',16)
l=legend('$\hat{W}_{c3,1}$','$\hat{W}_{c3,2}$','$\hat{W}_{c3,3}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPE1LinkWc3Nash3Player.pdf
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
l=legend('$\hat{W}_{a1,1}$','$\hat{W}_{a1,2}$','$\hat{W}_{a1,3}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPE1LinkWa1Nash3Player.pdf
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
l=legend('$\hat{W}_{a2,1}$','$\hat{W}_{a2,2}$','$\hat{W}_{a2,3}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPE1LinkWa2Nash3Player.pdf
clear p
close all

%% WaH2

p=plot(WaH3.time,WaH3.signals.values,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{W}_{a2}(t)$','FontSize',16,'Interpreter','latex')
title('Player 3 Actor Weights','FontSize',16)
set(gca,'FontSize',16)
l=legend('$\hat{W}_{a3,1}$','$\hat{W}_{a3,2}$','$\hat{W}_{a3,3}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPE1LinkWa3Nash3Player.pdf
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
l=legend('$x_1$','$x_2$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPE1LinkStateNash3Player.pdf
clear p
close all

%% Control

p=plot(u1.time,u1.signals.values,u2.time,u2.signals.values,u3.time,u3.signals.values,'LineWidth',2);
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
l=legend('Player 1','Player 2','Player 3');
set(l,'Location','southeast');
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPE1LinkControlNash3Player.pdf
clear p
close all

%% xTilde

p=plot(xTD.time,xTD.signals.values,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
axis([0 15 -1 1])
axis 'auto x'
xlabel('Time (s)','FontSize',16)
ylabel('$\dot{\tilde{x}}$','FontSize',16,'Interpreter','latex')
title('State Derivative Estimation Error','FontSize',16)
set(gca,'FontSize',16)
l=legend('$\dot{\tilde{x}}_1$','$\dot{\tilde{x}}_2$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPPE1LinkxTildeDotNash3Player.pdf
clear p
close all

close all
