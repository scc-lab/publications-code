%% Probing Signal Plots
n=4;
h=plot(t,z(:,1:4),'LineWidth',2);
title ('System States','FontSize',16);
xlabel ('Time (s)');
h1=legend('$q_1 \: (rad)$','$q_2 \: (rad)$','$\dot{q}_1 \: (rad/s)$','$\dot{q}_2 \: (rad/s)$');
set(h1,'interpreter','latex','Location','SouthEast')
axis([0 40 -Inf Inf])
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf JP1ADPTrackingPEStates.pdf
close all

h=plot(t,E,'LineWidth',2);
title ('Tracking Error','FontSize',16);
xlabel ('Time (s)');
h1=legend('$e_1 \: (rad)$','$e_2 \: (rad)$','$e_3 \: (rad/s)$','$e_4 \: (rad/s)$');
set(h1,'interpreter','latex','Location','SouthEast')
axis([0 40 -Inf Inf])
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf JP1ADPTrackingPEError.pdf
close all

plot(t,z(:,n+1:n+L),'LineWidth',2);
title ('Parameters of the critic NN','FontSize',16);
xlabel ('Time (s)');
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf JP1ADPTrackingPEWc.pdf
close all

plot(t,z(:,n+L+1:n+2*L),'LineWidth',2); 
title ('Parameters of the actor NN','FontSize',16);
xlabel ('Time (s)');
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf JP1ADPTrackingPEWa.pdf
close all


%% Comparison Plots

load 'JP1ADPTracking2LinkADP.mat'
load 'JP1ADPTracking2LinkGPOPS.mat'

%% Control Plot
h1=plot(output.solutionPlot.time,output.solutionPlot.control,'LineWidth',2);
grid on
hold on
h2 = plot(ADPsol.time,ADPsol.control,':','LineWidth',2);
set(cell2mat(get(cell2mat(get(h1,'Annotation')),'LegendInformation')),'IconDisplayStyle','off');
set(cell2mat(get(cell2mat(get(h2,'Annotation')),'LegendInformation')),'IconDisplayStyle','off');
a1=cell2mat(get(cell2mat(get(h1,'Annotation')),'LegendInformation'));
a2=cell2mat(get(cell2mat(get(h2,'Annotation')),'LegendInformation'));
set(a1(1),'IconDisplayStyle','on');
set(a2(1),'IconDisplayStyle','on');
n=100;
M_n = round(linspace(1,size(output.solutionPlot.control,1),n)) ;
m = {'*','s'}';
ms = {6,5}';
h3=plot(output.solutionPlot.time(M_n',:),output.solutionPlot.control(M_n',:),'s');
mc = get(h3,'color');
set(h3,{'marker'},m,{'markers'},ms,{'MarkerFaceColor'},mc);
n=20;
M_n = round(linspace(1,size(output.solutionPlot.control,1),n)) ;
h4=plot(ADPsol.time(M_n',:),ADPsol.control(M_n',:),'s');
set(h4,{'marker'},m,{'markers'},ms,{'MarkerFaceColor'},mc);
set(cell2mat(get(cell2mat(get(h4,'Annotation')),'LegendInformation')),'IconDisplayStyle','off');
axis([0 20 -Inf Inf])
set(gca,'FontSize',16)
xlabel('Time(s)')
title('\mu(\zeta)')
h1=legend('GPOPS','ADP','$\mu_1$','$\mu_2$');
set(h1,'interpreter','latex','Location','SouthEast')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf JP1ADPTrackingPolicy.pdf
close all
%% State Plots

h1=plot(output.solutionPlot.time,output.solutionPlot.state(:,1:4),'LineWidth',2);
grid on
hold on
h2 = plot(ADPsol.time,ADPsol.state,':','LineWidth',2);
set(cell2mat(get(cell2mat(get(h1,'Annotation')),'LegendInformation')),'IconDisplayStyle','off');
set(cell2mat(get(cell2mat(get(h2,'Annotation')),'LegendInformation')),'IconDisplayStyle','off');
a1=cell2mat(get(cell2mat(get(h1,'Annotation')),'LegendInformation'));
a2=cell2mat(get(cell2mat(get(h2,'Annotation')),'LegendInformation'));
set(a1(1),'IconDisplayStyle','on');
set(a2(1),'IconDisplayStyle','on');
n=100;
M_n = round(linspace(1,size(output.solutionPlot.state,1),n)) ;
m = {'o','^','s','*'}';
ms = {5,5,4,6}';
h3=plot(output.solutionPlot.time(M_n',:),output.solutionPlot.state(M_n',1:4),'s');
mc = get(h3,'color');
set(h3,{'marker'},m,{'markers'},ms,{'MarkerFaceColor'},mc);
n=20;
M_n = round(linspace(1,size(output.solutionPlot.state,1),n)) ;
h4=plot(ADPsol.time(M_n',:),ADPsol.state(M_n',:),'s');
set(h4,{'marker'},m,{'markers'},ms,{'MarkerFaceColor'},mc);
set(cell2mat(get(cell2mat(get(h4,'Annotation')),'LegendInformation')),'IconDisplayStyle','off');
axis([0 20 -Inf Inf])
set(gca,'FontSize',16)
xlabel('Time(s)')
title('Tracking Error')
h1=legend('GPOPS','ADP','$e_1$','$e_2$','$e_3$','$e_4$')
set(h1,'interpreter','latex','Location','SouthEast')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf JP1ADPTrackingError.pdf
close all
%% costate
h1=plot(output.solutionPlot.time,output.solutionPlot.costate,'LineWidth',1);
grid on
hold on
set(cell2mat(get(cell2mat(get(h1,'Annotation')),'LegendInformation')),'IconDisplayStyle','off');
n=100;
M_n = round(linspace(1,size(output.solutionPlot.costate,1),n)) ;
h3=plot(output.solutionPlot.time(M_n',:),output.solutionPlot.costate(M_n',:),'s');
m = {'o','^','s','*','x','d','v','+'}';
ms = {5,5,4,6,6,4,5,6}';
mc = get(h3,'color');
set(h3,{'marker'},m,{'markers'},ms,{'MarkerFaceColor'},mc);
axis([0 20 -Inf Inf])
set(gca,'FontSize',16)
xlabel('Time(s)')
title('Costate')
h4=legend('$\partial V^*/\partial e_1$','$\partial V^*/\partial e_2$','$\partial V^*/\partial e_3$','$\partial V^*/\partial e_4$',...
    '$\partial V^*/\partial x_{d1}$','$\partial V^*/\partial x_{d2}$','$\partial V^*/\partial x_{d3}$','$\partial V^*/\partial x_{d4}$');
set(h4,'FontSize',14,'Location','SouthEast','Interpreter','Latex');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf JP1ADPTrackingCostate.pdf
close all
%% Hamiltonian
plot(output.solution.time,output.solution.Hamiltonian,'LineWidth',2)
axis([0 20 -Inf Inf])
grid on
set(gca,'FontSize',16)
title('Hamiltonian From GPOPS');
xlabel('Time(s)')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf JP1ADPTrackingGPOPSHamilt.pdf
close all