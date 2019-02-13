function dae = OC2LinkGPOPSdae(soldae)

t = soldae.time;
Zeta = soldae.state;
mU = soldae.control;
p1=3.473;
p2=.196;
p3=.242;
fd1=5.3;
fd2=1.1;
fs1=8.45;
fs2=2.35;
Fd = diag([fd1,fd2]);
Fcl = zeros(length(t),8);
for i=1:length(t)
    zeta=Zeta(i,:);
    e = zeta(1:4);
    mu = mU(i,:);
    xd = zeta(4+1:end);
    h = [xd(3) xd(4) -4*xd(1) -9*xd(2)];
    x = e+xd;
    Minv = (1/(p2^2 - p1*p2 + p3^2*cos(x(2))^2))*...
       [-p2               p2+p3*cos(x(2))  ;...
         p2+p3*cos(x(2)) -p1-2*p3*cos(x(2))];
    Vm = [-p3*sin(x(2))*x(4) -p3*sin(x(2))*(x(3)+x(4));
           p3*sin(x(2))*x(3)  0                       ];
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
    ud = gplusd*(h'-fd);
    F = [f-h'+g*ud;h'];
    G = [g;zeros(size(g))];
    Fcl(i,:) = F'+mu*G';
end
path = [];
dae = [Fcl path];