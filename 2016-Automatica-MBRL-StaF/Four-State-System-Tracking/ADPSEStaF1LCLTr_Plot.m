%% WcH
WcH=getElement(ADPSEStaF1LCLTrData,'WcH');
plot(WcH.Values.Time(1:100:end,:),WcH.Values.Data(1:100:end,:),'LineWidth',2);
xlabel('Time (s)','FontSize',20)
ylabel('$\hat{W}_{c}(t)$','FontSize',20,'Interpreter','latex')
title('Value Function Weights','FontSize',20)
grid on
set(gca,'FontSize',20)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPSEStaF1LCLTrDataWc.pdf
close all
clear WcH
%% WaH
WaH=getElement(ADPSEStaF1LCLTrData,'WaH');
plot(WaH.Values.Time(1:100:end,:),WaH.Values.Data(1:100:end,:),'LineWidth',2);
xlabel('Time (s)','FontSize',20)
ylabel('$\hat{W}_{a}(t)$','FontSize',20,'Interpreter','latex')
title('Policy Weights','FontSize',20)
grid on
set(gca,'FontSize',20)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPSEStaF1LCLTrDataWa.pdf
close all
clear WaH
%% Theta
ThetaH=getElement(ADPSEStaF1LCLTrData,'SIDThetaHat');
Tt=ThetaH.Values.Time(1:100:end,:);
TD=reshape(ThetaH.Values.Data(:,:,1:100:end),6,1001);
clear ThetaH
p=plot(Tt,TD,'LineWidth',2);
grid on
mrk1={'s','v','o','*','^','+'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',20)
ylabel('$\hat{\theta}(t)$','FontSize',18,'Interpreter','latex')
title('Drift Dynamics NN Weights','FontSize',20)
set(gca,'FontSize',20)
hold on
grid on
p1=plot(Tt(1:100:end),-1*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(1),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
p1=plot(Tt(1:100:end),1*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(2),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
p1=plot(Tt(1:100:end),0*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(3),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
p1=plot(Tt(1:100:end),-0.5*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(4),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
p1=plot(Tt(1:100:end),0*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(5),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
p1=plot(Tt(1:100:end),-0.5*ones(size(Tt(1:100:end),1)),'--','LineWidth',2);
set(p1,'Color',get(p(6),'Color'),'LineWidth',2)
% hasbehavior(p1, 'legend', false);
l=legend('$\hat{\theta}_1$','$\hat{\theta}_2$','$\hat{\theta}_3$','$\hat{\theta}_4$','$\hat{\theta}_5$','$\hat{\theta}_6$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPSEStaF1LCLTrDataTheta.pdf
clear p Tt TD
close all

%% State
x=getElement(ADPSEStaF1LCLTrData,'e');
p=plot(x.Values.Time(1:100:end,:),x.Values.Data(1:100:end,:),'LineWidth',2);
xlabel('Time (s)','FontSize',20)
ylabel('$x(t)$','FontSize',20,'Interpreter','latex')
title('State Trajectory','FontSize',20)
grid on
set(gca,'FontSize',20)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPSEStaF1LCLTrDataState.pdf
clear p x
close all

%% Control
u=getElement(ADPSEStaF1LCLTrData,'u');
p=plot(u.Values.Time(1:100:end,:),u.Values.Data(1:100:end,:),'LineWidth',2);
xlabel('Time (s)','FontSize',20)
ylabel('$u(t)$','FontSize',20,'Interpreter','latex')
title('Control Trajectory','FontSize',20)
grid on
set(gca,'FontSize',20)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPSEStaF1LCLTrDataControl.pdf
clear p u
close all

%% Tracking Error

e=getElement(ADPSEStaF1LCLTrData,'e');
p=plot(e.Values.Time(1:100:end,:),e.Values.Data(1:100:end,:),'LineWidth',2);
xlabel('Time (s)','FontSize',20)
ylabel('$e(t)$','FontSize',20,'Interpreter','latex')
title('Tracking Error','FontSize',20)
grid on
set(gca,'FontSize',20)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPSEStaF1LCLTrDataError.pdf
clear p e
close all

%% ADP regressor matrix minimum eigenvalue

ADPeig=getElement(ADPSEStaF1LCLTrData,'ADP_eig');
p=plot(ADPeig.Values.Time(1:100:end,:),ADPeig.Values.Data(1:100:end,:),'LineWidth',2);
xlabel('Time (s)','FontSize',20)
title('Minimum Eigenvalue of BE Extrapolation Matrix','FontSize',20)
grid on
set(gca,'FontSize',20)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPSEStaF1LCLTrDataADPEig.pdf
clear p ADPeig
close all

%% SysID regressor matrix minimum singular value

CLeig=getElement(ADPSEStaF1LCLTrData,'SIDRank');
p=plot(CLeig.Values.Time(1:100:end,:),CLeig.Values.Data(1:100:end,:),'LineWidth',2);
xlabel('Time (s)','FontSize',20)
title('Minimum Singular Value of CL History Stack','FontSize',20)
grid on
set(gca,'FontSize',20)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPSEStaF1LCLTrDataCLEig.pdf
clear p CLeig
close all