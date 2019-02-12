%% State Phase
xd11=[0;0;0];
xd31=[-0.5;-0.5;0];
xd21=[0.5;-0.5;0];
xd23=xd21-xd31;
xd40=[0.5;0.5;0];
xd50=[-0.5;0.5;0];
x0=getElement(ADPSEMREMKNNData,'x0');
mrk1={'+', 'o', '*' , 'x' , 'square' , 'diamond' , 'v' , '^' , '>' , '<' , 'pentagram' , 'hexagram' };
col1={'b','k','r','g','m'};
for i=1:5
    c=num2str(i);
    if i<=3
        d='1';
    else
        d='0';
    end
    cc=num2str(i+1);
    eval(['x=getElement(ADPSEMREMKNNData,''x' c ''')']);
    eval(['plot3(x.Values.Time(1:100:end),x.Values.Data(1:100:end,1),x.Values.Data(1:100:end,2),col1{' c '},''LineWidth'',1)']);
    hold on;
    grid on
    eval(['p=plot3(x0.Values.Time(1:100:end),x0.Values.Data(1:100:end,1)+xd' c d '(1),x0.Values.Data(1:100:end,2)+xd' c d '(2),col1{' c '},''LineWidth'',0.75)']);
    eval(['p.Marker=mrk1{' cc '}']);
    p.MarkerFaceColor=p.Color;
    p.MarkerSize=6;
    f_nummarkers3(p,20,0);
    for j=1:size(p,1)
        hasbehavior(p(j), 'legend', false);
    end
end

xlabel('Time (s)','FontSize',16)
ylabel('$x_{i1}(t)$','FontSize',16,'Interpreter','latex')
zlabel('$x_{i2}(t)$','FontSize',16,'Interpreter','latex')
title('State Trajectories','FontSize',16)
l=legend('$x_1$','$x_0$','$x_2$','$x_0+x_{d20}$','$x_3$','$x_0+x_{d30}$','$x_4$','$x_0+x_{d40}$','$x_5$','$x_0+x_{d50}$');
set(l,'Interpreter','latex','Location','northeast');
grid on
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 5]);
set(gcf, 'PaperPosition', [0 0 8 5]);
print -dpdf ADPSEMREMKNNDataState.pdf
clear p x1 x2 x3
close all

%% State Phase
x0=getElement(ADPSEMREMKNNData,'x0');
for i=1:5
    c=num2str(i)
    eval(['x=getElement(ADPSEMREMKNNData,''x' c ''')']);
    eval(['plot(x.Values.Data(end-32000:end,1),x.Values.Data(end-32000:end,2),col1{' c '},''LineWidth'',2)']);
    hold on;
end
plot(x0.Values.Data(end-32000:end,1),x0.Values.Data(end-32000:end,2),'g--','LineWidth',2);
xlabel('$x(t)$','FontSize',20,'Interpreter','latex')
ylabel('$y(t)$','FontSize',20,'Interpreter','latex')
title('State Trajectories','FontSize',20)
grid on
set(gca,'FontSize',20)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPSEMREMKNNDataState2D.pdf
clear p x1 x2 x3
close all

%% Control
for i=1:5
    c=num2str(i);
    eval(['u=getElement(ADPSEMREMKNNData,''u' c ''')']);
    plot(u.Values.Time(1:100:end),u.Values.Data(1:100:end,:),'LineWidth',2);
    xlabel('Time (s)','FontSize',20)
    eval(['ylabel(''$u_' c '(t)$'',''FontSize'',20,''Interpreter'',''latex'')']);
    title('Control Trajectory','FontSize',20)
    grid on
    set(gca,'FontSize',20)
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 5]);
    set(gcf, 'PaperPosition', [0 0 6 5]);
    eval(['print -dpdf ADPSEMREMKNNDataControl' c '.pdf']);
    clear u
    close all
end

%% Tracking Error
for i=1:5
    c=num2str(i);
    eval(['e=getElement(ADPSEMREMKNNData,''e' c ''')']);
    plot(e.Values.Time(1:100:end),e.Values.Data(1:100:end,:),'LineWidth',2);
    xlabel('Time (s)','FontSize',20)
    eval(['ylabel(''$e_' c '(t)$'',''FontSize'',20,''Interpreter'',''latex'')']);
    title('Tracking Error','FontSize',20)
    grid on
    set(gca,'FontSize',20)
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 5]);
    set(gcf, 'PaperPosition', [0 0 6 5]);
    eval(['print -dpdf ADPSEMREMKNNDataError' c '.pdf']);
    clear p e
    close all
end

%% WcH
for i=1:5
    c=num2str(i);
    eval(['WcH=getElement(ADPSEMREMKNNData,''WcH' c ''')']);
    plot(WcH.Values.Time(1:100:end),WcH.Values.Data(1:100:end,:),'LineWidth',2);
    xlabel('Time (s)','FontSize',20)
    eval(['ylabel(''$\hat{W}_{c' c '}\:(t)$'',''FontSize'',20,''Interpreter'',''latex'')']);
    title('Value Function Weights','FontSize',20)
    grid on
    set(gca,'FontSize',20)
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 5]);
    set(gcf, 'PaperPosition', [0 0 6 5]);
    eval(['print -dpdf ADPSEMREMKNNDataWcH' c '.pdf'])
end
close all
%% WaH
for i=1:5
    c=num2str(i);
    eval(['WaH=getElement(ADPSEMREMKNNData,''WaH' c ''')']);
    plot(WaH.Values.Time(1:100:end),WaH.Values.Data(1:100:end,:),'LineWidth',2);
    xlabel('Time (s)','FontSize',20)
    eval(['ylabel(''$\hat{W}_{a' c '}\:(t)$'',''FontSize'',20,''Interpreter'',''latex'')']);
    title('Policy Weights','FontSize',20)
    grid on
    set(gca,'FontSize',20)
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 5]);
    set(gcf, 'PaperPosition', [0 0 6 5]);
    eval(['print -dpdf ADPSEMREMKNNDataWaH' c '.pdf'])
end
close all

% %% ADP regressor matrix minimum eigenvalue
% 
% ADPeig=getElement(ADPSEMREMKNNData,'ADP_eig');
% p=plot(ADPeig.Values.Time(1:100:end),ADPeig.Values.Data,'LineWidth',2);
% xlabel('Time (s)','FontSize',20)
% title('Minimum Eigenvalue of BE Extrapolation Matrix','FontSize',20)
% grid on
% set(gca,'FontSize',20)
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [6 5]);
% set(gcf, 'PaperPosition', [0 0 6 5]);
% print -dpdf ADPSEMREMKNNDataADPEig.pdf
% clear p ADPeig
% close all
% 
% %% SysID regressor matrix minimum singular value
% 
% CLeig=getElement(ADPSEMREMKNNData,'SIDRank');
% p=plot(CLeig.Values.Time(1:100:end),CLeig.Values.Data,'LineWidth',2);
% xlabel('Time (s)','FontSize',20)
% title('Minimum Singular Value of CL History Stack','FontSize',20)
% grid on
% set(gca,'FontSize',20)
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [6 5]);
% set(gcf, 'PaperPosition', [0 0 6 5]);
% print -dpdf ADPSEMREMKNNDataCLEig.pdf
% clear p CLeig
% close all