function [c,ceq] = OC_2linkrobot_const(x)

OrderOfCollocation=45;
dimx=4;
dimu=2;
F=zeros(OrderOfCollocation,dimx);
ceq2 = zeros(OrderOfCollocation,dimx);

D=OCGenD(OrderOfCollocation);
[a,w]=OClgwt(OrderOfCollocation,-1,1);
a = sort(a);

X=vec2mat(x(1:((OrderOfCollocation+1)*dimx)),OrderOfCollocation+1,dimx); X=X';
U=vec2mat(x(((OrderOfCollocation+1)*dimx+1):end),OrderOfCollocation,dimu); U=U';

p1=3.473;
p2=.196;
p3=.242;
fd1=5.3;
fd2=1.1;
fs1=8.45;
fs2=2.35;

tempxinf = 0;
for i=1:OrderOfCollocation
    
    tempx=X(i+1,:);
    tempu=U(i,:);
    M  = [ p1+2*p3*cos(tempx(2))   , p2+p3*cos(tempx(2));
           p2+p3*cos(tempx(2))     , p2              ];
    Vm = [ -p3*sin(tempx(2))*tempx(4)  , -p3*sin(tempx(2))*(tempx(3)+tempx(4));
            p3*sin(tempx(2))*tempx(3)  , 0                                 ];
    Fd = diag([fd1,fd2]);
    Fs = [fs1*tanh(tempx(3));fs2*tanh(tempx(4))];
    
    F(i,:) = [tempx(3) tempx(4) [inv(M)*((-Vm - Fd)*[tempx(3); tempx(4)]+[tempu(1);tempu(2)]-Fs)]'];
    
    ceq2(i,:) = D(i,:)*X - (2/(1-a(i))^2)*F(i,:);
    tempxinf = tempxinf + w(i)*(2/(1-a(i))^2)*F(i,:);    
%     ceq2(i,:) = D(i,:)*X - (1/(1-a(i)))*F(i,:);
%     tempxinf = tempxinf + w(i)*(1/(1-a(i)))*F(i,:);
%     ceq2(i,:) = D(i,:)*X - (2/(1-a(i)))*F(i,:);
%     tempxinf = tempxinf + w(i)*(2/(1-a(i)))*F(i,:);

end

x0=[1 1 0 0]';
xinf = tempxinf + x0';

ceq1=[X(1,:)'-x0; xinf'];

ceq2=reshape(ceq2',OrderOfCollocation*dimx,1);
ceq=[ceq1;ceq2];

c=[];

