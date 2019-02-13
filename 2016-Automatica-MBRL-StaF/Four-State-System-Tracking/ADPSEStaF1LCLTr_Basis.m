function [sig,sig_p,L]=ADPSEStaF1LCLTr_Basis(Z,s,n)
%#codegen
    d1=s.*[0;0;1/sqrt(3);1];
    d2=s.*[0;0;1/sqrt(3);-1];
    d3=s.*[0;0;-2/sqrt(3);0];
    d4=s.*[0;0;1/sqrt(10);1/sqrt(6)];
    d5=s.*[0;0;-2*sqrt(2/5);0];
    d6=[0;0;0;0];
    e=Z(1:n,1);
    xd=Z(n+1:end,1);
    sig=[exp(Z'*(Z+d1))-1;
         exp(Z'*(Z+d2))-1;
         exp(Z'*(Z+d3))-1;
         exp(Z'*(Z+d4))-1;
         exp(Z'*(Z+d5))-1];
%          exp(Z'*(Z+d6))-1];
    sig_p=[(2*Z+d1)'*exp(Z'*(Z+d1));
           (2*Z+d2)'*exp(Z'*(Z+d2));
           (2*Z+d3)'*exp(Z'*(Z+d3));
           (2*Z+d4)'*exp(Z'*(Z+d4));
           (2*Z+d5)'*exp(Z'*(Z+d5))];
%            (2*Z+d6)'*exp(Z'*(Z+d6))];
%     sig=[exp(Z'*(Z+d1))-exp(xd'*(xd+d1(n+1:end,1)));
%          exp(Z'*(Z+d2))-exp(xd'*(xd+d2(n+1:end,1)));
%          exp(Z'*(Z+d3))-exp(xd'*(xd+d3(n+1:end,1)));
%          exp(Z'*(Z+d4))-exp(xd'*(xd+d4(n+1:end,1)));
%          exp(Z'*(Z+d5))-exp(xd'*(xd+d5(n+1:end,1)))];
%     sig_p=[(2*Z+d1)'*exp(Z'*(Z+d1))-[zeros(1,n) (2*xd+d1(n+1:end,1))'*exp(xd'*(xd+d1(n+1:end,1)))];
%            (2*Z+d2)'*exp(Z'*(Z+d2))-[zeros(1,n) (2*xd+d2(n+1:end,1))'*exp(xd'*(xd+d2(n+1:end,1)))];
%            (2*Z+d3)'*exp(Z'*(Z+d3))-[zeros(1,n) (2*xd+d3(n+1:end,1))'*exp(xd'*(xd+d3(n+1:end,1)))];
%            (2*Z+d4)'*exp(Z'*(Z+d4))-[zeros(1,n) (2*xd+d4(n+1:end,1))'*exp(xd'*(xd+d4(n+1:end,1)))];
%            (2*Z+d5)'*exp(Z'*(Z+d5))-[zeros(1,n) (2*xd+d5(n+1:end,1))'*exp(xd'*(xd+d5(n+1:end,1)))]];
   L=size(sig_p,1);
end

             
             
