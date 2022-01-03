function [ W ] = rk4method( x_init, f, h, T)
%RK4 Method



W = x_init;

for i=1:floor(T/h)-1
    V1 = f(W(:,i));
    V2 = f(W(:,i) + 1/2*h*V1);
    V3 = f(W(:,i) + 1/2*h*V2);
    V4 = f(W(:,i) + h*V3);
    
    W = [W,W(:,i) + h/6*(V1+2*V2+2*V3+V4)];
end


end

