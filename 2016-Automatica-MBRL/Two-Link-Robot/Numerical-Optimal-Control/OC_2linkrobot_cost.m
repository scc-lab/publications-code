function J = OC_2linkrobot_cost(x)

OrderOfCollocation = 45;
[a,w]=OClgwt(OrderOfCollocation,-1,1);
a = sort(a);
temp=0;

dimx=4;
dimu=2;

% Q = 10*eye(dimx);.
Q=[10 0 0 0; 0 10 0 0; 0 0 1 0; 0 0 0 1];
R = eye(dimu);

X=vec2mat(x(1:((OrderOfCollocation+1)*dimx)),OrderOfCollocation+1,dimx); X=X';
U=vec2mat(x(((OrderOfCollocation+1)*dimx+1):end),OrderOfCollocation,dimu); U=U';

for k = 1:OrderOfCollocation
    tempx=X(k+1,:);
    tempu=U(k,:);
    L=tempx*Q*tempx' + tempu*R*tempu';
     
    temp=temp+w(k)*(2/(1-a(k))^2)*L;
%     temp=temp+w(k)*(1/(1-a(k)))*L;
%     temp=temp+w(k)*(2/(1-a(k)))*L;

end
J=temp;