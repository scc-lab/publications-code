%% WcH

WcH=getElement(ADPTLinearCL,'WcH');
Tt=WcH.Values.Time;
TDtemp=WcH.Values.Data;
TD=TDtemp(:,[1 2 5 3 4 6:end]);
clear WcH
p=plot(Tt,TD,'LineWidth',2);
mrk1={'s','v','o','none','none','none','none','none','none','none','none','none','none','none','none'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{W}_c(t)$','FontSize',14,'Interpreter','latex')
title('Value function NN Weights','FontSize',16)
set(gca,'FontSize',16)
hold on
p1=plot(Tt(1:100:end),4.4331*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(1),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
p1=plot(Tt(1:100:end),1.3526*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(2),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
p1=plot(Tt(1:100:end),2.9036*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(3),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
l=legend('$\hat{W}_{c1}$','$\hat{W}_{c2}$','$\hat{W}_{c5}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPTLinearCLWcH.pdf
clear p Tt TD TDtemp
close all

%% WaH
WaH=getElement(ADPTLinearCL,'WaH');
Tt=WaH.Values.Time;
TDtemp=WaH.Values.Data;
TD=TDtemp(:,[1 2 5 3 4 6:end]);
clear WaH
p=plot(Tt,TD,'LineWidth',2);
mrk1={'s','v','o','none','none','none','none','none','none','none','none','none','none','none','none'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{W}_a(t)$','FontSize',14,'Interpreter','latex')
title('Policy NN Weights','FontSize',16)
set(gca,'FontSize',16)
hold on
p1=plot(Tt(1:100:end),4.4331*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(1),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
p1=plot(Tt(1:100:end),1.3526*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(2),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
p1=plot(Tt(1:100:end),2.9036*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(3),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
l=legend('$\hat{W}_{a1}$','$\hat{W}_{a2}$','$\hat{W}_{a5}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPTLinearCLWaH.pdf
clear p Tt TD TDtemp
close all

%% Theta
ThetaH=getElement(ADPTLinearCL,'SIDThetaHat');
Tt=ThetaH.Values.Time;
TD=reshape(ThetaH.Values.Data,4,10001);
clear ThetaH
p=plot(Tt,TD,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\hat{\theta}(t)$','FontSize',14,'Interpreter','latex')
title('Drift Dynamics NN Weights','FontSize',16)
set(gca,'FontSize',16)
hold on
p1=plot(Tt(1:100:end),(-1)*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(1),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
p1=plot(Tt(1:100:end),1*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(2),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
p1=plot(Tt(1:100:end),(-0.5)*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(3),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
p1=plot(Tt(1:100:end),(-0.5)*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(4),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
l=legend('$\hat{\theta}(1,1)$','$\hat{\theta}(2,1)$','$\hat{\theta}(1,2)$','$\hat{\theta}(2,2)$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPTLinearCLTheta.pdf
clear p Tt TD
close all

%% State
x=getElement(ADPTLinearCL,'x');
p=plot(x.Values.Time,x.Values.Data,'LineWidth',2);
xlabel('Time (s)','FontSize',16)
ylabel('$x(t)$','FontSize',16,'Interpreter','latex')
title('State Trajectory','FontSize',16)
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPTLinearCLState.pdf
clear p x
close all

%% Control
u=getElement(ADPTLinearCL,'u');
p=plot(u.Values.Time,u.Values.Data,'LineWidth',2);
xlabel('Time (s)','FontSize',16)
ylabel('$u(t)$','FontSize',16,'Interpreter','latex')
title('Control Trajectory','FontSize',16)
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPTLinearCLControl.pdf
clear p u
close all

%% Tracking Error

e=getElement(ADPTLinearCL,'e');
p=plot(e.Values.Time,e.Values.Data,'LineWidth',2);
xlabel('Time (s)','FontSize',16)
ylabel('$e(t)$','FontSize',16,'Interpreter','latex')
title('Tracking Error','FontSize',16)
grid on
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPTLinearCLError.pdf
clear p e
close all

%% ADP regressor matrix minimum eigenvalue

ADPeig=getElement(ADPTLinearCL,'ADP_eig');
p=plot(ADPeig.Values.Time,ADPeig.Values.Data,'LineWidth',2);
xlabel('Time (s)','FontSize',16)
title('Minimum Eigenvalue of BE Extrapolation Matrix','FontSize',16)
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPTLinearCLADPEig.pdf
clear p ADPeig
close all

%% SysID regressor matrix minimum singular value

CLeig=getElement(ADPTLinearCL,'SIDRank');
p=plot(CLeig.Values.Time,CLeig.Values.Data,'LineWidth',2);
xlabel('Time (s)','FontSize',16)
title('Minimum Singular Value of CL History Stack','FontSize',16)
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPTLinearCLCLEig.pdf
clear p CLeig
close all