function [sig_p,L]=ADPSELSCLTr_Basis(Z,adp)
%#codegen
e1=Z(1);e2=Z(2);xd1=Z(3);xd2=Z(4);
if adp==1
    sig=[e1^2;e1*e2;e1*xd1;e1*xd2;e2^2;e2*xd1;e2*xd2;xd1^2;xd2^2;xd1*xd2];
sig_p =[[ 2*e1,    0,  0,  0];
        [   e2,   e1,  0,  0];
        [  xd1,    0, e1,  0];
        [  xd2,    0,  0, e1];
        [    0, 2*e2,  0,  0];
        [    0,  xd1, e2,  0];
        [    0,  xd2,  0, e2]];
L=size(sig_p,1);
else
sig_p =[[ 2*e1,    0,     0,     0];
        [   e2,   e1,     0,     0];
        [  xd1,    0,    e1,     0];
        [  xd2,    0,     0,    e1];
        [    0, 2*e2,     0,     0];
        [    0,  xd1,    e2,     0];
        [    0,  xd2,     0,    e2];
        [    0,    0, 2*xd1,     0];
        [    0,    0,     0, 2*xd2];
        [    0,    0,   xd2,   xd1]];
L=size(sig_p,1);
end

             
             
