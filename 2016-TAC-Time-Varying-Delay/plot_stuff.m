figure(1)
clf
subplot(3,1,1);
hold on;
plot(error_case1.time,error_case1.signals.values(:,1),'--r','LineWidth',2);
plot(error_case1.time,error_case1.signals.values(:,2),'-k','LineWidth',2);
ylabel('Error','FontSize',12,'FontName','Times');
leg1 = legend('e_1','e_2');
set(leg1,'FontSize',7);
%axis([0 60 -2e-3 2.5e-3]);
grid on
subplot(3,1,2);
hold on;
plot(u1_case1.time,u1_case1.signals.values(:,1),'--r','LineWidth',2);
plot(u1_case1.time,u1_case1.signals.values(:,2),'-k','LineWidth',2);
ylabel('Control','FontSize',12,'FontName','Times');
leg1 = legend('u_1','u_2');
set(leg1,'FontSize',7);
%axis([0 60 -2e-3 2.5e-3]);
grid on
subplot(3,1,3);
hold on;
plot(tau_case1.time,1000*tau_case1.signals.values(:,1),'--r','LineWidth',2);
plot(tau_case1.time,1000*tau_case1.signals.values(:,2),'-k','LineWidth',2);
xlabel('Time (s)','FontSize',12,'FontName','Times');
ylabel('Delay (ms)','FontSize',12,'FontName','Times');
leg1 = legend('tau_s','tau_i');
set(leg1,'FontSize',7);
%axis([0 60 -2e-3 2.5e-3]);
grid on
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
%print -dpdf case5.pdf;