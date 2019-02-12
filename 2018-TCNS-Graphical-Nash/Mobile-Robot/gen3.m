function [ mu3,scrG3,scrF3,calG3,calF3,Gsig3,F3 ] = gen3( calE3,calE1,Wa3H,R3 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
e1 = calE1(1:3);
x1 = calE1(4:end);
x3 = calE3(7:9);
x0 = x1-e1;
[ ~, f10 ] = gen0( x0 );
g1=geng(x1);
g3=geng(x3);
g31 = eye(2);
f31 = [0;0];
% u31 = f31+g31*u1;
scrLg3=[eye(2) -g31;zeros(2) eye(2)];
scrLg3inv=[eye(2) eye(2);zeros(2) eye(2)];
scrLg33=[eye(2) eye(2)];
scrLg31=[zeros(2) eye(2)];
F3 = [f31;f10];
scrF3=g3*scrLg33*F3-g1*scrLg31*F3;
scrG3=g3*scrLg33-g1*scrLg31;
calF3=g3*scrLg33*F3;
calG3=g3*scrLg33;
[~,sig3GRe3,sig3GRx3,~,~]=ADPSEMREMKNN_Basis3(calE3);         
Gsig3=(scrG3*[eye(2);zeros(2)])'*sig3GRe3'+(calG3*[eye(2);zeros(2)])'*sig3GRx3';
mu3=-0.5*(R3\Gsig3)*Wa3H;
end

