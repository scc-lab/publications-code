%% WcH
for i=1:N
    c=num2str(i);
    LL=cumsum(L);
    WcHi=WcH(LL(i)-L(i)+1:LL(i),:);
    plot(T,WcHi,'LineWidth',2);
    xlabel('Time (s)','FontSize',16)
    eval(['ylabel(''$\hat{W}_{c' c '}\:(t)$'',''FontSize'',16,''Interpreter'',''latex'')']);
    title('Critic weights','FontSize',16)
    axis([0 30 -10 10]);
    axis 'auto y'
    set(gca,'FontSize',16)
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 5]);
    set(gcf, 'PaperPosition', [0 0 6 5]);
    % print -dpdf ADPPEFree2LinkWc.pdf
    eval(['print -dpdf ADPCLNNWcH' c '.pdf'])
end
close all

%% WaH

for i=1:N
    c=num2str(i);
    LL=cumsum(L);
    WaHi=WaH(LL(i)-L(i)+1:LL(i),:);
    plot(T,WaHi,'LineWidth',2);
    xlabel('Time (s)','FontSize',16)
    eval(['ylabel(''$\hat{W}_{a' c '}(t)$'',''FontSize'',16,''Interpreter'',''latex'')']);
    title('Actor weights','FontSize',16)
    axis([0 30 -10 10]);
    axis 'auto y'
    set(gca,'FontSize',16)
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 5]);
    set(gcf, 'PaperPosition', [0 0 6 5]);
    % print -dpdf ADPPEFree2LinkWc.pdf
    eval(['print -dpdf ADPCLNNWaH' c '.pdf'])
end
close all
%% State

p=plot(T,X,'LineWidth',2);
mrk1={'s','v','o','*','^','+'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time (s)','FontSize',16)
ylabel('$x(t)$','FontSize',16,'Interpreter','latex')
title('State trajectories','FontSize',16)
set(gca,'FontSize',16)
hold on
lcell=cell(1,N);
for i=1:N
p1=plot(T(1:100:end),X0(i,1:100:end),'--','LineWidth',2);
set(p1,'Color',get(p(i),'Color'),'LineWidth',2)
eval(['lcell{' num2str(i) '}=''$x_' num2str(i) '$'';']);
end
l=legend(lcell{:});
axis([0 30 -10 10]);
    axis 'auto y'
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLNNX.pdf
clear p
close all
%% Error

plot(T(2:end),E,'LineWidth',2);
set(gca,'FontSize',16);
title('Error trajectories');
for i=1:N
eval(['lcell{' num2str(i) '}=''$e_' num2str(i) '$'';']);
end
h=legend(lcell{:});
set(h,'Interpreter','latex','Location','northeast')
xlabel('Time (s)','FontSize',16)
ylabel('$e(t)$','Interpreter','latex','FontSize',16);
axis([0 30 -10 10]);
    axis 'auto y'
grid on
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLNNE.pdf
close all
%% Control

plot(T(2:end),U,'LineWidth',2);
set(gca,'FontSize',16);
title('Control trajectories');
for i=1:N
eval(['lcell{' num2str(i) '}=''$u_' num2str(i) '$'';']);
end
h=legend(lcell{:});
set(h,'Interpreter','latex','Location','northeast')
xlabel('Time (s)','FontSize',16)
ylabel('$u(t)$','Interpreter','latex','FontSize',16);
grid on
axis([0 30 -10 1])
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLNNU.pdf
close all

%% Virtual Control

plot(T(2:end),MU,'LineWidth',2);
set(gca,'FontSize',16);
title('Relative control error rajectories');
for i=1:N
eval(['lcell{' num2str(i) '}=''$\mu_' num2str(i) '$'';']);
end
h=legend(lcell{:});
set(h,'Interpreter','latex','Location','northeast')
xlabel('Time (s)','FontSize',16)
ylabel('$\mu(t)$','Interpreter','latex','FontSize',16);
axis([0 30 -10 10]);
    axis 'auto y'
grid on
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
print -dpdf ADPCLNNMU.pdf
close all

%% ThetaH
for i=1:N
    c=num2str(i);
    PP=cumsum(P);
    ThetaHi=ThetaH(PP(i)-P(i)+1:PP(i),:);
    p=plot(T,ThetaHi,'LineWidth',2);
    mrk1={'s','v','o','*','^','+'};
    mrk=(mrk1(1,1:size(p,1)))';
    set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',10);
    f_nummarkers(p,10,1);
    for j=1:size(p,1)
    hasbehavior(p(j), 'legend', false);
    end
    xlabel('Time (s)','FontSize',16)
    eval(['ylabel(''$\hat{\theta}_{' c '}\:(t)$'',''FontSize'',16,''Interpreter'',''latex'')']);
    title('Drift parameters','FontSize',16)
    set(gca,'FontSize',16)
    hold on
    if i==1, Th=[0;1]; elseif i==2, Th=[0,0.5]; elseif i==3, Th=[0.1,1]; elseif i==4, Th=[0.2,1]; else Th=[0.5,1]; end
    p1=plot(T(1:100:end),Th(1)*ones(size(T(1:100:end),2)),'--','LineWidth',2);
    set(p1,'Color',get(p(1),'Color'),'LineWidth',2)
    % hasbehavior(p1, 'legend', false);
    p1=plot(T(1:100:end),Th(2)*ones(size(T(1:100:end),2)),'--','LineWidth',2);
    set(p1,'Color',get(p(2),'Color'),'LineWidth',2)
    eval(['l=legend(''$\hat{\theta}_{' c '1}$'',''$\hat{\theta}_{' c '2}$'');']);
    set(l,'Interpreter','latex','Location','northeast')
%     axis([0 tf -0.5 1.5])
    axis([0 30 -10 10]);
    axis 'auto y'
    set(gca,'FontSize',16)
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 5]);
    set(gcf, 'PaperPosition', [0 0 6 5]);
    eval(['print -dpdf ADPCLNNThetaH' c '.pdf'])
    close all
end

%% Total Cost 

plot(T(2:end),Lintegral,'linewidth',2)
xlabel('Time (s)','FontSize',16);
ylabel('$J(t)$','Interpreter','latex','FontSize',16);
for i=1:N
eval(['lcell{' num2str(i) '}=''$J_' num2str(i) '$'';']);
end
h=legend(lcell{:});
set(h,'Interpreter','latex','Location','northeast')
title('Accumulated Cost','FontSize',16)
axis([0 30 -10 10]);
axis 'auto y'
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
eval(['print -dpdf ADPCLNNJ.pdf'])
close all
