% Rushikesh Kamalapurkar
% ADP Tracking

%% Initialization
clear all
clc

%% Problem Definition
test = 0;
tf = 60;
n = 4; % Number of States
R = 1; 
Q = [10 0  0  0;
     0  10 0  0;
     0  0  2  0;
     0  0  0  2];

%% Control Gains
pecutoff = 0.5;
basis = 1; % 1 = Known basis OR 2 = Symmetric Sigmoids.
if basis == 1
    L = 23;
else
    L = 10;
end
Vc = -100+200*rand(2*n+1,L);
etaa1 = 5;
etaa2 = 0.001;
lambda = 0.001;
etac = 1.25;

v = 0.005;

%% Initial Conditions stack
% x0 = [0.8;0.6;0;0];
x0 = [1.8;1.6;0;0];
WcH = 10*ones(L,1);
WaH = 6*ones(L,1);
if test == 1
    WcH = [83.3664360398476;2.37005961853738;27.0535257208229;2.78151919205775;-2.82560210947881;0.201121140905491;14.1349718360624;29.8140596906650;18.8766797804977;4.11019727718757;3.47631801536084;6.69761592266286;9.71065522881670;15.5795959278823;4.97416255918046;12.4247304839404;11.3106634637451;3.29418471671678;1.19713547876155;-1.99945807222058;4.55196662573689;-0.468991113555371;0.565476966194629];
    WaH = WcH;
    etaa1 = 0;
    etaa2 = 0;
    etac = 0;
    lambda = 0;
    pecutoff = 0;
end
Gamma = reshape(2000*eye(L),L*L,1);
z0=[x0;WcH;WaH;Gamma];
 
%% Integration
options = odeset('OutputFcn',@odeplot,'OutputSel',[1 2 3 4],'RelTol',...
    1e-3, 'AbsTol', 1e-3);
figure(5);
[t,z] = ode23(@(t,z)ADP_Tracking_2Link_dynamics(t,z,basis,tf,n,R,Q,...
    Vc,L,etaa1,etaa2,etac,lambda,v,pecutoff),...
    [0 tf],z0,options);
% fs = 500;
% t = 0:1/fs:tf;
% [z]= ode1(@(t,z)ADP_Tracking_2Link_dynamics(t,z,basis,tf,n,R,Q,...
%     Vc,L,etaa1,etaa2,etac,v,pecutoff),...
%     t,z0);

%% For Plotting
Xd = zeros(length(t),n);
E = zeros(length(t),n);
OMEGA = zeros(length(t),L);
DELTA = zeros(length(t),1);
MU = zeros(length(t),2);
U = zeros(length(t),2);
Cost = zeros(size(t));
for i=1:length(t)
    x = z(i,1:n)';
    xd = [(1/2)*cos(2*t(i)); (1/3)*cos(3*t(i)); -sin(2*t(i)); -sin(3*t(i))];
%     h = [xd(3); xd(4); -4*xd(1); -9*xd(2)];
%     xd = [(1/4)*sin(2*t(i)); (1/3)*cos(3*t(i)); (1/2)*cos(2*t(i)); -sin(3*t(i))];
    h = [xd(3); xd(4); -4*xd(1); -9*xd(2)];
    e = x-xd;
    E(i,:)=e';
    Xd(i,:)= xd'; 
    Z = [e;xd];
    WcH = z(i,n+1:n+L)';
    WaH = z(i,n+L+1:n+2*L)';
    Gamma = reshape(z(i,n+2*L+1:n+2*L+L*L),L,L);
    p1 = 3.473;
    p2 = .196;
    p3 = .242;
    fd1 = 5.3;
    fd2 = 1.1;
    fs1 = 8.45;
    fs2 = 2.35;
    M = [p1+2*p3*cos(x(2)) p2+p3*cos(x(2));
         p2+p3*cos(x(2))   p2             ];
    Minv = (1/(p2^2 - p1*p2 + p3^2*cos(x(2))^2))*...
           [-p2               p2+p3*cos(x(2))  ;...
             p2+p3*cos(x(2)) -p1+2*p3*cos(x(2))];
    Vm = [-p3*sin(x(2))*x(4) -p3*sin(x(2))*(x(3)+x(4));
           p3*sin(x(2))*x(3)  0                       ];
    Fd = diag([fd1,fd2]);
    Fs = [fs1*tanh(x(3));fs2*tanh(x(4))];
    f = [x(3); x(4); Minv*((-Vm - Fd)*[x(3); x(4)]-Fs)];
    g = [0 0; 0 0; Minv]; 
    Md = [p1+2*p3*cos(xd(2)) p2+p3*cos(xd(2));
          p2+p3*cos(xd(2))   p2              ];
    Vmd = [-p3*sin(xd(2))*xd(4) -p3*sin(xd(2))*(xd(3)+xd(4));
            p3*sin(xd(2))*xd(3)  0                          ];
    Fsd = [fs1*tanh(xd(3));fs2*tanh(xd(4))];
    fd = [xd(3); xd(4); Md\((-Vmd - Fd)*[xd(3); xd(4)]-Fsd)];
    gplusd = [0 0; 0 0; Md']';
    ud = gplusd*(h-fd);
    F = [f-h+g*ud;h];
    G = [g;zeros(size(g))];
    if basis == 1
        phi = [(1/2)*Z(1)^2;...
               (1/2)*Z(2)^2;...
               (1/2)*Z(3)^2;...
               (1/2)*Z(4)^2;...
               Z(1)*Z(3);...
               Z(1)*Z(4);...
               Z(2)*Z(3);...
               Z(2)*Z(4);...
               (1/2)*Z(1)^2*Z(2)^2;...
               (1/2)*Z(1)^2*Z(5)^2;...
               (1/2)*Z(1)^2*Z(6)^2;...
               (1/2)*Z(1)^2*Z(7)^2;...
               (1/2)*Z(1)^2*Z(8)^2;...
               (1/2)*Z(2)^2*Z(5)^2;...
               (1/2)*Z(2)^2*Z(6)^2;...
               (1/2)*Z(2)^2*Z(7)^2;...
               (1/2)*Z(2)^2*Z(8)^2;...
               (1/2)*Z(3)^2*Z(4)^2;...
               (1/2)*Z(3)^2*Z(5)^2;...
               (1/2)*Z(3)^2*Z(6)^2;...
               (1/2)*Z(3)^2*Z(7)^2;...
               (1/2)*Z(3)^2*Z(8)^2;...
               (1/2)*Z(4)^2*Z(5)^2;...
               (1/2)*Z(4)^2*Z(6)^2;...
               (1/2)*Z(4)^2*Z(7)^2;...
               (1/2)*Z(4)^2*Z(8)^2];

        phi_p = [Z(1)   0      0      0      0      0      0      0      ;...        
                 0      Z(2)   0      0      0      0      0      0      ;...   
%                  0      0      Z(3)   0      0      0      0      0      ;...
%                  0      0      0      Z(4)   0      0      0      0      ;...
                 Z(3)   0      Z(1)   0      0      0      0      0      ;...
                 Z(4)   0      0      Z(1)   0      0      0      0      ;...
                 0      Z(3)   Z(2)   0      0      0      0      0      ;...        
                 0      Z(4)   0      Z(2)   0      0      0      0      ;...
                 Z(1)*Z(2)^2 Z(2)*Z(1)^2 0           0           0           0           0           0           ;...
                 Z(1)*Z(5)^2 0           0           0           Z(5)*Z(1)^2 0           0           0           ;...  
                 Z(1)*Z(6)^2 0           0           0           0           Z(6)*Z(1)^2 0           0           ;...
                 Z(1)*Z(7)^2 0           0           0           0           0           Z(7)*Z(1)^2 0           ;...
                 Z(1)*Z(8)^2 0           0           0           0           0           0           Z(8)*Z(1)^2 ;...
                 0           Z(2)*Z(5)^2 0           0           Z(5)*Z(2)^2 0           0           0           ;...  
                 0           Z(2)*Z(6)^2 0           0           0           Z(6)*Z(2)^2 0           0           ;... 
                 0           Z(2)*Z(7)^2 0           0           0           0           Z(7)*Z(2)^2 0           ;...
                 0           Z(2)*Z(8)^2 0           0           0           0           0           Z(8)*Z(2)^2 ;...
%                  0           0           Z(3)*Z(4)^2 Z(4)*Z(3)^2 0           0           0           0           ;...
                 0           0           Z(3)*Z(5)^2 0           Z(5)*Z(3)^2 0           0           0           ;...
                 0           0           Z(3)*Z(6)^2 0           0           Z(6)*Z(3)^2 0           0           ;...
                 0           0           Z(3)*Z(7)^2 0           0           0           Z(7)*Z(3)^2 0           ;...
                 0           0           Z(3)*Z(8)^2 0           0           0           0           Z(8)*Z(3)^2 ;...
                 0           0           0           Z(4)*Z(5)^2 Z(5)*Z(4)^2 0           0           0           ;...
                 0           0           0           Z(4)*Z(6)^2 0           Z(6)*Z(4)^2 0           0           ;...
                 0           0           0           Z(4)*Z(7)^2 0           0           Z(7)*Z(4)^2 0           ;...
                 0           0           0           Z(4)*Z(8)^2 0           0           0           Z(8)*Z(4)^2 ];
    else
        phi=(1-exp(-2*(Vc'*x)))./(1+exp(-2*(Vc'*x)));
        phi_p = diag(4*exp(-2*Vc'*[1;Z])./(1+exp(-2*Vc'*[1;Z])).^2)*Vc(2:end,:)';
    end
    mu=-0.5*(R\G')*phi_p'*WaH;
    pe = 0.85*[tanh(2*t(i))*(20*sin(sqrt(230)*pi*t(i))*cos(sqrt(20)*pi*t(i)) + 6*sin(18*exp(2)*t(i))+ 20*cos(40*t(i))*cos(21*t(i)));...
         0];
    if t<=tf*pecutoff
        munew = mu + pe;
    else
        munew = mu;
    end
    u = munew+ud;
    omega = phi_p*(F+G*mu);
    OMEGA(i,:) = omega;
    r = e'*Q*e + mu'*R*mu;
    Cost(i) = r;
    delta = WcH'*omega + r;
    DELTA(i,:) = delta;
    U(i,:) = u;
    MU(i,:) = munew;
    xD = f+g*u;
end
cost = trapz(t,Cost);
ADPsol.cost = cost;
ADPsol.Hamiltonian = DELTA;
ADPsol.control = MU;
ADPsol.policy = U;
ADPsol.state = E;
ADPsol.time = t;
%% State Plots
fig1 = figure (6);
subplot(2,3,1)
plot(t,z(:,1:4));
title ('System States');
xlabel ('Time (s)');
subplot(2,3,2)
plot(t,z(:,n+1:n+L));
title ('Parameters of the critic NN');
xlabel ('Time (s)');
legend ('W_{c1}','W_{c2}', 'W_{c3}');
subplot(2,3,3)
plot(t,z(:,n+L+1:n+2*L)); 
title ('Parameters of the actor NN');
xlabel ('Time (s)');
legend ('W_{a1}','W_{a2}', 'W_{a3}');
subplot(2,3,4)
plot(t,MU);
grid on
title ('Munew');
xlabel ('Time (s)');
subplot(2,3,5)
plot(t,E); 
title ('Tracking Error');
grid on
xlabel ('Time (s)');
subplot(2,3,6)
plot(t,DELTA); 
title ('Regressor');
xlabel ('Time(s)');
if test == 1
    clearvars -except ADPsol
    figure(2)
    hold on
    plot(ADPsol.time,ADPsol.control,'-');
    grid on
    title ('mu');
    xlabel ('Time (s)');
    axis([0,30,-Inf,Inf]);
    figure(3)
    hold on
    plot(ADPsol.time,ADPsol.policy,'-');
    grid on
    title ('u');
    xlabel ('Time (s)');
    axis([0,30,-Inf,Inf]);
    figure(1)
    hold on
    plot(ADPsol.time,ADPsol.state,'-'); 
    title ('Tracking Error');
    grid on
    xlabel ('Time (s)');
    axis([0,30,-Inf,Inf]);
    figure(4)
    hold on
    plot(ADPsol.time,ADPsol.Hamiltonian,'-'); 
    title ('Bellman Error');
    grid on
    xlabel ('Time (s)');
    axis([0,30,-Inf,Inf]);
end

    