clear all
clc

dimx=4;
dimu=2;
OrderOfCollocation=45;

Z0x1 = 0*ones(OrderOfCollocation+1,1); 
Z0x2 = 1*ones(OrderOfCollocation+1,1);
Z0x3 = 0*ones(OrderOfCollocation+1,1);
Z0x4 = 1*ones(OrderOfCollocation+1,1);
Z0u1 = 2*ones(OrderOfCollocation,1);
Z0u2 = 1*ones(OrderOfCollocation,1);
Z0=[Z0x1;Z0x2;Z0x3;Z0x4;Z0u1;Z0u2];

o=optimset('MaxFunEvals',1000000,'TolCon',1e-8,...
    'Display','iter','TolFun',1e-5,'Algorithm','interior-point');
[sol,~,e,~,lambda] = fmincon(@OC_2linkrobot_cost,Z0,[],[],[],[],[],[],@OC_2linkrobot_const,o);

X=vec2mat(sol(1:((OrderOfCollocation+1)*dimx)),OrderOfCollocation+1,dimx); 
X=X';
U=vec2mat(sol(((OrderOfCollocation+1)*dimx+1):end),OrderOfCollocation,dimu);
U=U';

[a,w]=OClgwt(OrderOfCollocation,-1,1);
a = sort(a);

dtau = 0.001;
Tau = -1:dtau:(1-dtau);
Control = zeros(length(Tau),dimu);
for tau = -1:dtau:(1-dtau);
    for i = 1:OrderOfCollocation
        L = 1;
        for j = 1:OrderOfCollocation
            if j~=i
                L = L*((tau-a(j))/(a(i)-a(j)));
            end
        end
        Control(nearest(tau*(1/dtau)+(1/dtau)+1),:) = ...
            Control(nearest(tau*(1/dtau)+(1/dtau)+1),:) + U(i,:)*L;
    end
end


td = (ones(size(a))+a)./(ones(size(a))-a);
tc = (ones(size(Tau))+Tau)./(ones(size(Tau))-Tau);
% td = log(2./(ones(size(a))-a));
% tc = log(2./(ones(size(Tau))-Tau));
% td = log(4./(ones(size(a))-a).^2);
% tc = log(4./(ones(size(Tau))-Tau).^2);


plot(td,U,'*')
hold on
plot(tc,Control)
grid on
figure;
plot(td,X(2:end,:),'*')
grid on;
hold on;