function [sigth,p]=ADPSE1LCLTr_SIDBasis(x,cl,bias,Ynn)
%#codegen
if cl==1
    if bias==0
        sigth=tanh(Ynn'*x);
    else
        sigth=[1;tanh(Ynn'*x)];
    end
elseif cl==2
    if bias==0
        sigth=[x(1);x(2);x(1)*x(2);x(1)^2;x(1)^2;x(1)^2*x(2);x(1)*x(2)^2;x(1)^3;x(2)^3];
    else
        sigth=[1;x(1);x(2);x(1)*x(2);x(1)^2;x(1)^2;x(1)^2*x(2);x(1)*x(2)^2;x(1)^3;x(2)^3];
    end
else
    sigth=[x(1);x(2);x(2)*(1-(cos(2*x(1))+2)^2)];
end
p=size(sigth,1);