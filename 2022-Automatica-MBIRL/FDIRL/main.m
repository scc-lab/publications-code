clear all
close all
evalc("linearPreSim");
sim("linearSim.slx");

%% Tracking Error
x=getElement(logsout,'e');
T=x.Values.Time;
X=squeeze(x.Values.Data);
p=plot(T,X,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$e(t)$','FontSize',16,'Interpreter','latex')
set(gca,'FontSize',16)
l=legend('$e_{(1)}$','$e_{(2)}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 3]);
set(gcf, 'PaperPosition', [0 0 6 3]);
print -dpdf LinNoNoiseError.pdf
clear p
close all

%% Parameter Estimation Error
x=getElement(logsout,'theta_tilde');
T=x.Values.Time;
X=x.Values.Data;
p=plot(T,X,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\tilde{\theta}(t)$','FontSize',16,'Interpreter','latex')
set(gca,'FontSize',16)
l=legend('$\tilde{\theta}_{(1)}$','$\tilde{\theta}_{(2)}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 3]);
set(gcf, 'PaperPosition', [0 0 6 3]);
print -dpdf LinNoNoiseThetaTilde.pdf
clear p
close all

%% State Estimation Error
x=getElement(logsout,'x_tilde');
T=x.Values.Time;
X=x.Values.Data;
p=plot(T,X,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\tilde{x}(t)$','FontSize',16,'Interpreter','latex')
set(gca,'FontSize',16)
l=legend('$\tilde{x}_{(1)}$','$\tilde{x}_{(2)}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 3]);
set(gcf, 'PaperPosition', [0 0 6 3]);
print -dpdf LinNoNoisexTilde.pdf
clear p
close all

%% Feedback Estimation Error
x=getElement(logsout,'W_u_tilde');
T=x.Values.Time;
X=squeeze(x.Values.Data);
p=plot(T,X,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\tilde{W}_u(t)$','FontSize',16,'Interpreter','latex')
set(gca,'FontSize',16)
l=legend('$\tilde{W}_{u(1)}$','$\tilde{W}_{u(2)}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 3]);
set(gcf, 'PaperPosition', [0 0 6 3]);
print -dpdf LinNoNoiseWuTilde.pdf
clear p
close all

%% IRL Estimation Error
x=getElement(logsout,'W_tilde');
T=x.Values.Time;
X=squeeze(x.Values.Data);
p=plot(T,X,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\tilde{W}(t)$','FontSize',16,'Interpreter','latex')
set(gca,'FontSize',16)
l=legend('$\tilde{W}_{(1)}$','$\tilde{W}_{(2)}$','$\tilde{W}_{(3)}$','$\tilde{W}_{(4)}$','$\tilde{W}_{(5)}$');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 3]);
set(gcf, 'PaperPosition', [0 0 6 3]);
print -dpdf LinNoNoiseWTilde.pdf
clear p l x T X mrk mrk1 i ans
close all

save(['Data - ' strrep(datestr(now), ':', '-') '.mat'])
%% Simulation Without feedback estimation
IRL.feedback_driven=0;
sim("linearSim.slx");

%% IRL Estimation Error
x=getElement(logsout,'W_tilde');
T=x.Values.Time;
X=squeeze(x.Values.Data);
p=plot(T,X,'LineWidth',2);
mrk1={'s','v','o','*','^'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$\tilde{W}(t)$','FontSize',16,'Interpreter','latex')
set(gca,'FontSize',16)
l=legend('$\tilde{W}_{(1)}$','$\tilde{W}_{(2)}$','$\tilde{W}_{(3)}$','$\tilde{W}_{(4)}$','$\tilde{W}_{(5)}$');
set(l,'Interpreter','latex','Location','southwest');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 3]);
set(gcf, 'PaperPosition', [0 0 6 3]);
print -dpdf LinNoNoiseWTildeNoWu.pdf
clear p
close all
