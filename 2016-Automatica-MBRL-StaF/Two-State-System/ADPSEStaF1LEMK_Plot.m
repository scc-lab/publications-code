Refsize = 3;
%% WcH
WcHT=getElement(ADPCLStaF,'WcH');
TT=WcHT.Values.Time(1:100:end,:);
WcH=WcHT.Values.Data(1:100:end,:);
p=plot(TT,WcH,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',20)
ylabel('$\hat{W}_{c}(t)$','FontSize',20,'Interpreter','latex')
title('Value Function Weights','FontSize',20)
set(gca,'FontSize',20)
l=legend('$\hat{W}_{c,1}$','$\hat{W}_{c,2}$','$\hat{W}_{c,3}$');
grid on
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLStaf1LWcH.pdf
clear p
close all

%% WaH

WaHT=getElement(ADPCLStaF,'WaH');
TT=WaHT.Values.Time(1:100:end,:);
WaH=WaHT.Values.Data(1:100:end,:);
p=plot(TT,WaH,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',20)
ylabel('$\hat{W}_{a}(t)$','FontSize',20,'Interpreter','latex')
title('Policy Weights','FontSize',20)
set(gca,'FontSize',20)
l=legend('$\hat{W}_{a,1}$','$\hat{W}_{a,2}$','$\hat{W}_{a,3}$');
grid on
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLStaf1LWaH.pdf
clear p
close all

%% State

XT=getElement(ADPCLStaF,'x');
TT=XT.Values.Time(1:100:end,:);
X=XT.Values.Data(1:100:end,:);
p=plot(TT,X,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',20)
ylabel('$x(t)$','FontSize',20,'Interpreter','latex')
title('State Trajectory','FontSize',20)
set(gca,'FontSize',20)
l=legend('$x_1$','$x_2$');
grid on
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLStaf1LState.pdf
clear p
close all

%% Control

UT=getElement(ADPCLStaF,'u');
TT=UT.Values.Time(1:100:end,:);
U=UT.Values.Data(1:100:end,:);
p=plot(TT,U,'LineWidth',2);
xlabel('Time (s)','FontSize',20)
ylabel('$u\left(t\right)$','FontSize',20,'Interpreter','latex')
title('Control Trajectory','FontSize',20)
grid on
% l=legend('Player 1','Player 2');
% set(l,'Location','southeast');
set(gca,'FontSize',20)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLStaf1LControl.pdf
clear p
close all

%% Value function estimation

VVT=getElement(ADPCLStaF,'VError');
TT=VVT.Values.Time(1:100:end,:);
VV=VVT.Values.Data(1:100:end,:);
p=plot(TT,VV,'LineWidth',2);
xlabel('Time (s)','FontSize',20)
h=ylabel('$V^{*}\left(x(t)\right)-\hat{V}\left(x(t),\hat{W}_c(t)\right)$','FontSize',18,'Interpreter','latex');
set(h,'unit','character')
title('Value Function Estimation Error','FontSize',20)

grid on
% l=legend('Player 1','Player 2');
% set(l,'Location','southeast');
set(gca,'FontSize',20)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLStaf1LVStarmV.pdf
clear p
close all

%% Control Estimation Error

VVT=getElement(ADPCLStaF,'uError');
TT=VVT.Values.Time(1:100:end,:);
VV=VVT.Values.Data(1:100:end,:);
p=plot(TT,VV,'LineWidth',2);
xlabel('Time (s)','FontSize',20)
h=ylabel('$u^{*}\left(x(t)\right)-\hat{u}\left(x(t),\hat{W}_a(t)\right)$','FontSize',18,'Interpreter','latex');
set(h,'unit','character')
title('Control Estimation Error','FontSize',20)
grid on
% l=legend('Player 1','Player 2');
% set(l,'Location','southeast');
set(gca,'FontSize',20)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLStaf1LuError.pdf
clear p
close all
%% Policy estimation

GradVST=getElement(ADPCLStaF,'uStar');
GradVHT=getElement(ADPCLStaF,'u');
TT=GradVST.Values.Time(1:100:end,:);
GradVS=GradVST.Values.Data(1:100:end,:);
GradVH=GradVHT.Values.Data(1:100:end,:);
p=plot(TT(1:end),[GradVS(1:end,:) GradVH(1:end,:)],'LineWidth',1.5);
mrk1={'s','v','o','^','*'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
axis([0 5 -2 0.1])
l=legend('$u^*(x(t))$','$\hat{u}(x(t),\hat{W}_a(t))$');
grid on
set(l,'Interpreter','latex','Location','southeast');
xlabel('Time (s)','FontSize',20)
% h=ylabel('$\nabla V \left(x(t)\right)$','FontSize',20,'Interpreter','latex');
% set(h,'unit','character');
title('Optimal Control Estimation','FontSize',20)
grid on
% l=legend('Player 1','Player 2');
% set(l,'Location','southeast');
set(gca,'FontSize',20)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLStaf1LOptCont.pdf
clear p
% close all

%% Minimum Singular Value of CL Matrix
CBART=getElement(ADPCLStaF,'cbar');
TT=CBART.Values.Time(1:100:end,:);
CBAR=CBART.Values.Data(1:100:end,:);
p=plot(TT,CBAR,'LineWidth',2);
xlabel('Time (s)','FontSize',20)
ylabel('$\underline{c}(t)$','FontSize',20,'Interpreter','latex')
title('Minimum Singular Value of CL Matrix','FontSize',20)
grid on
% l=legend('Player 1','Player 2');
% set(l,'Location','southeast');
set(gca,'FontSize',20)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLStaf1Lcbar.pdf
clear p
close all