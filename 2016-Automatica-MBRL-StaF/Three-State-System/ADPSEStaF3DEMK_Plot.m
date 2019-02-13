Refsize = 3;
%% WcH
WcHT=getElement(ADPCLStaF,'WcH');
TT=WcHT.Values.Time;
WcH=WcHT.Values.Data;
p=plot(TT,WcH,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{W}_{c}(t)$','FontSize',16,'Interpreter','latex')
title('Value Function Weights','FontSize',16)
set(gca,'FontSize',16)
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
TT=WaHT.Values.Time;
WaH=WaHT.Values.Data;
p=plot(TT,WaH,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{W}_{a}(t)$','FontSize',16,'Interpreter','latex')
title('Policy Weights','FontSize',16)
set(gca,'FontSize',16)
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
TT=XT.Values.Time;
X=XT.Values.Data;
p=plot(TT,X,'LineWidth',2);
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
TT=UT.Values.Time;
U=UT.Values.Data;
p=plot(TT,U,'LineWidth',2);
xlabel('Time (s)','FontSize',16)
ylabel('$u\left(t\right)$','FontSize',16,'Interpreter','latex')
title('Control Trajectory','FontSize',16)
grid on
% l=legend('Player 1','Player 2');
% set(l,'Location','southeast');
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLStaf1LControl.pdf
clear p
close all

%% Value function estimation

VVT=getElement(ADPCLStaF,'VError');
TT=VVT.Values.Time;
VV=VVT.Values.Data;
p=plot(TT,VV,'LineWidth',2);
xlabel('Time (s)','FontSize',16)
h=ylabel('$V^{*}\left(x(t)\right)-\hat{V}\left(x(t),\hat{W}_c(t)\right)$','FontSize',14,'Interpreter','latex');
set(h,'unit','character')
title('Value Function Estimation Error','FontSize',16)
grid on
% l=legend('Player 1','Player 2');
% set(l,'Location','southeast');
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLStaf1LVStarmV.pdf
clear p
close all

%% Control Estimation Error

VVT=getElement(ADPCLStaF,'uError');
TT=VVT.Values.Time;
VV=VVT.Values.Data;
p=plot(TT,VV,'LineWidth',2);
xlabel('Time (s)','FontSize',16)
h=ylabel('$u^{*}\left(x(t)\right)-\hat{u}\left(x(t),\hat{W}_a(t)\right)$','FontSize',14,'Interpreter','latex');
set(h,'unit','character')
title('Control Estimation Error','FontSize',16)
grid on
% l=legend('Player 1','Player 2');
% set(l,'Location','southeast');
set(gca,'FontSize',16)
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
TT=GradVST.Values.Time;
GradVS=GradVST.Values.Data;
GradVH=GradVHT.Values.Data;
p=plot(TT(1:end),[GradVS(1:end,:) GradVH(1:end,:)],'LineWidth',1.5);
mrk1={'s','v','o','^','*'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
axis([0 10 -2 0.1])
l=legend('$u^*(x(t))$','$\hat{u}(x(t),\hat{W}_a(t))$');
grid on
set(l,'Interpreter','latex','Location','southeast');
xlabel('Time (s)','FontSize',16)
% h=ylabel('$\nabla V \left(x(t)\right)$','FontSize',16,'Interpreter','latex');
% set(h,'unit','character');
title('Optimal Control Estimation','FontSize',16)
grid on
% l=legend('Player 1','Player 2');
% set(l,'Location','southeast');
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLStaf1LOptCont.pdf
clear p
close all

%% Minimum Singular Value of CL Matrix
CBART=getElement(ADPCLStaF,'cbar');
TT=CBART.Values.Time;
CBAR=CBART.Values.Data;
p=plot(TT,CBAR,'LineWidth',2);
xlabel('Time (s)','FontSize',16)
ylabel('$\underline{c}(t)$','FontSize',16,'Interpreter','latex')
title('Minimum Singular Value of CL Matrix','FontSize',16)
grid on
% l=legend('Player 1','Player 2');
% set(l,'Location','southeast');
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLStaf1Lcbar.pdf
clear p
close all