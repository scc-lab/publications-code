function el_sim
clc

saveData = false;
paperFigs = false;

p1 = 3.473; p2 = 0.196; p3 = 0.242;
fd1 = 5.3; fd2 = 1.1;
theta = [p1; p2; p3; fd1; fd2];
x0 = [zeros(4,1); 10*rand(5,1); zeros(2,1); zeros(2*5,1)];
alpha = 1;
k1 = 1; k2 = 3;
Gamma = 0.3*eye(5);
delT = 0.1;
tspan = [0,20];
initTime = 1; % must be >= delT
eigThresh = 0.01;
N = 20;

% Run once (short tspan) to get history for DDE and prefill history stacks
sol = dde23(@sys,delT,x0(1:9),linspace(tspan(1),initTime,3000)',[],0,[],[],theta,alpha,k1,k2,Gamma,eigThresh);
t = sol.x';
x = sol.y';
dx = zeros(9,length(t));
tau = zeros(2,length(t));
Y4 = zeros(2,5,length(t));
for ii = 1:length(t)
    [dx(:,ii),tau(:,ii),Y4(:,:,ii)] = sys(t(ii),x(ii,:)',[],0,[],[],theta,alpha,k1,k2,Gamma,eigThresh);
end
dx = dx';
tau = tau';
Y4 = reshape(Y4,2*5,length(t))';

% Calculate history of script states
scriptU = zeros(length(t),2);
dU = zeros(length(t),2);
scriptY4 = zeros(length(t),2*5);
dY4 = zeros(length(t),2*5);
for ii = 1:length(t)
    intTime = linspace(max(t(ii)-delT,0),t(ii),100)';
    if intTime(end)-intTime(1) > 0
        thisTau = interp1(t,tau,intTime);
        thisY4 = interp1(t,Y4,intTime);
        dU(ii,:) = thisTau(end,:) - thisTau(1,:);
        dY4(ii,:) = thisY4(end,:) - thisY4(1,:);
        scriptU(ii,:) = trapz(intTime,thisTau);
        scriptY4(ii,:) = trapz(intTime,thisY4);
    end
end
sol.y = [x,scriptU,scriptY4]';
sol.yp = [dx,dU,dY4]';

% Fill history stackes
Ystack = zeros(2*N,5);
Ustack = zeros(2*N,1);
t = linspace(2*delT,initTime,N)';
tDel = max(t-delT,0);
x = deval(sol,t)';
xDel = deval(sol,tDel)';
for ii = 1:N
    [~,~,~,scriptU,scriptY] = sys(t(ii),x(ii,:)',xDel(ii,:)',delT,zeros(2*N,5),zeros(2*N,1),theta,alpha,k1,k2,Gamma,eigThresh);
    Ustack(2*ii-1:2*ii) = scriptU;
    Ystack(2*ii-1:2*ii,:) = scriptY;
end

% Repeatedly run until sim end time reached.
% When event triggered (indicating new data to add to history stack), stop,
% add data, and restart sim.
done = false;
options = ddeset('Events',@events);
while ~done
    newSol = dde23(@sys,delT,sol,[sol.x(end),tspan(2)],options,delT,Ystack,Ustack,theta,alpha,k1,k2,Gamma,eigThresh);
    sol.x = [sol.x, newSol.x];
    sol.y = [sol.y, newSol.y];
    sol.yp = [sol.yp, newSol.yp];
    sol.x(end)
    
    x = deval(sol,[sol.x(end)-delT; sol.x(end)])';
    [~,~,~,scriptU,scriptY] = sys(sol.x(end),x(end,:)',x(1,:)',delT,Ystack,Ustack,theta,alpha,k1,k2,Gamma,eigThresh);
    if any(all(reshape(Ustack,[2,1,N]) == 0,1))
        ind = find(all(reshape(Ustack,[2,1,N]) == 0,1),1);
        Ustack(2*ind-1:2*ind) = scriptU;
        Ystack(2*ind-1:2*ind,:) = scriptY;
    else
        value = zeros(N,1);
        sum = Ystack'*Ystack;
        thisSum = scriptY'*scriptY;
        currEig = min(eig(sum));
        tot = sum+thisSum;
        for ii = 1:N
            value(ii) = min(eig(tot - Ystack(2*ii-1:2*ii,:)'*Ystack(2*ii-1:2*ii,:))) - currEig - eigThresh;
        end
        [~,ind] = max(value)
        Ustack(2*ind-1:2*ind) = scriptU;
        Ystack(2*ind-1:2*ind,:) = scriptY;
    end
    
    if sol.x(end) >= tspan(2)
        done = true;
    end
end
tICL = linspace(tspan(1),tspan(2),200)';
xICL = deval(sol,tICL)';

q = xICL(:,1:2);
qDot = xICL(:,3:4);
thetaHat = xICL(:,5:9);

[qd,qdDot,qdDotDot] = qDes(tICL');
qd = qd';
qdDot = qdDot';
qdDotDot = qdDotDot';

figure(1)
subplot(2,2,1)
plot(tICL,q*180/pi,tICL,qd*180/pi);
title('ICL States')
ylabel('ICL States')
ylim([-200,200])

subplot(2,2,2)
plot(tICL,bsxfun(@minus,theta',thetaHat));
title('ICL Est')
ylabel('ICL Est')
ylim([-10,10])

[tGrad,xGrad] = ode45(@sys,tspan,x0(1:9),[],x0,0,Ystack,Ustack,theta,alpha,k1,0,Gamma,eigThresh);

q = xGrad(:,1:2);
qDot = xGrad(:,3:4);
thetaHat = xGrad(:,5:9);

[qd,qdDot,qdDotDot] = qDes(tGrad');
qd = qd';
qdDot = qdDot';
qdDotDot = qdDotDot';

subplot(2,2,3)
plot(tGrad,q*180/pi,tGrad,qd*180/pi);
title('Basic state')
ylabel('Basic state')
ylim([-200,200])

subplot(2,2,4)
plot(tGrad,bsxfun(@minus,theta',thetaHat));
title('Basic Est')
ylabel('Basic Est')
ylim([-10,10])

if saveData
    saveStruct = struct('tICL',tICL,'xICL',xICL,'tGrad',tGrad,'xGrad',xGrad);
    save('EL_sim_data','-struct','saveStruct');
end

end

function [dx,tau,Y4,scriptU,scriptY] = sys(t,x,xDel,delT,Ystack,Ustack,theta,alpha,k1,k2,Gamma,~)

q = x(1:2,:);
qDot = x(3:4,:);
thetaHat = x(5:9,:);
if delT ~= 0
    scriptU = x(10:11,:);
    scriptY4 = reshape(x(12:21,:),2,5);
end

[qd,qdDot,~] = qDes(t);

if delT ~= 0
    tDel = max(t-delT,0);
    qDel = xDel(1:2,:);
    qDotDel = xDel(3:4,:);
    thetaHatDel = xDel(5:9,:);
    [qdDel,qdDotDel,~] = qDes(tDel);
end

e = qd - q;
eDot = qdDot - qDot;
r = eDot + alpha*e;

if delT ~= 0
    eDel = qdDel - qDel;
    eDotDel = qdDotDel - qDotDel;
    rDel = eDotDel + alpha*eDel;
end

Y2 = calcY2(t,q,qDot,alpha);
tau = Y2*thetaHat + e + k1*r;

[Mnt,Vmnt,Fdnt,Mdotnt] = terms(q,qDot,qDot,qDot);
Y4 = -Mdotnt + Vmnt + Fdnt;

if delT ~= 0
    Y2Del = calcY2(tDel,qDel,qDotDel,alpha);
    tauDel = Y2Del*thetaHatDel + eDel + k1*rDel;
    
    [MntDel,VmntDel,FdntDel,MdotntDel] = terms(qDel,qDotDel,qDotDel,qDotDel);
    Y3 = Mnt - MntDel;
    Y4Del = -MdotntDel + VmntDel + FdntDel;
    scriptY = Y3+scriptY4;
end

[M,Vm,Fd] = dyn(q,qDot,theta);
if delT ~= 0
    thetaHatDot = Gamma*Y2'*r + k2*Gamma*Ystack'*(Ustack - Ystack*thetaHat);
    scriptY4dot = Y4-Y4Del;
    dx = [qDot; M\(tau - Vm*qDot - Fd*qDot); thetaHatDot; tau-tauDel; scriptY4dot(:)];
else
    thetaHatDot = Gamma*Y2'*r;
    dx = [qDot; M\(tau - Vm*qDot - Fd*qDot); thetaHatDot];
end

end

function [value,isterminal,direction] = events(t,x,xDel,delT,Ystack,Ustack,theta,alpha,k1,k2,Gamma,eigThresh)

N = size(Ustack,1)/2;
if any(all(reshape(Ustack,[2,1,N]) == 0,1))
    value = 0;
    isterminal = 1;
    direction = 0;
else
    value = zeros(N,1);
    isterminal = ones(N,1);
    direction = ones(N,1);
    [~,~,~,~,scriptY] = sys(t,x,xDel,delT,Ystack,Ustack,theta,alpha,k1,k2,Gamma,eigThresh);
    sum = Ystack'*Ystack;
    thisSum = scriptY'*scriptY;
    currEig = min(eig(sum));
    tot = sum+thisSum;
    for ii = 1:N
        value(ii) = min(eig(tot - Ystack(2*ii-1:2*ii,:)'*Ystack(2*ii-1:2*ii,:))) - currEig - eigThresh;
    end
    
end

end

function [M,Vm,Fd,Mdot] = terms(q,qDot,Mmult,Vmult)

c2 = cos(q(2));
s2 = sin(q(2));

M = [Mmult(1), Mmult(2), 2*c2*Mmult(1)+c2*Mmult(2), 0, 0;
     0, Mmult(1)+Mmult(2), c2*Mmult(1), 0, 0];
Vm = [0, 0, -s2*qDot(2)*Vmult(1)-s2*(qDot(1)+qDot(2))*Vmult(2), 0, 0;
      0, 0, s2*qDot(1)*Vmult(1), 0, 0];
Fd = [0, 0, 0, qDot(1), 0;
      0, 0, 0, 0, qDot(2)];
Mdot = [0, 0, -2*s2*qDot(1)*qDot(2)-s2*qDot(2)^2, 0, 0;
        0, 0, -s2*qDot(1)*qDot(2), 0, 0];

end
    
function Y2 = calcY2(t,q,qDot,alpha)

[qd,qdDot,qdDotDot] = qDes(t);
e = qd - q;
eDot = qdDot - qDot;

[Mnt,Vmnt,Fdnt] = terms(q,qDot,qdDotDot+alpha*eDot,qdDot+alpha*e);
Y2 = Mnt+Vmnt+Fdnt;

end

function [M,Vm,Fd,Mdot] = dyn(q,qDot,theta)
p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
fd1 = theta(4);
fd2 = theta(5);

c2 = cos(q(2));
s2 = sin(q(2));
c2Dot = -1*sin(q(2))*qDot(2);

M = [p1+2*p3*c2, p2+p3*c2;
     p2+p3*c2, p2];
Vm = [-p3*s2*qDot(2), -p3*s2*(qDot(1)+qDot(2));
      p3*s2*qDot(1), 0];
Fd = [fd1, 0;
      0, fd2];
Mdot = [2*p3*c2Dot, p3*c2Dot;
        p3*c2Dot, 0];
end

function [qd,qdDot,qdDotDot] = qDes(t)

qd = [(1+10*exp(-2*t)).*sin(t);
      (1+10*exp(-t)).*cos(3*t)]*0.3;

qdDot = [cos(t).*(10*exp(-2*t) + 1) - 20*exp(-2*t).*sin(t);
         -3*sin(3*t).*(10*exp(-t) + 1) - 10*cos(3*t).*exp(-t)]*0.3;

qdDotDot = [40*exp(-2*t).*sin(t) - 40*exp(-2*t).*cos(t) - sin(t).*(10*exp(-2*t) + 1);
            10*cos(3*t).*exp(-t) - 9*cos(3*t).*(10*exp(-t) + 1) + 60*sin(3*t).*exp(-t)]*0.3;
end

