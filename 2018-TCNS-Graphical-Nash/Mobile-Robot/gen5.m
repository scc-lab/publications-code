function [ mu5,scrG5,scrF5,calG5,calF5,Gsig5,F5 ] = gen5( calE5,calE4,Wa5H,R5,xd40,xd45 )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
e4 = calE4(1:3);
x4 = calE4(end-2:end);
x5 = calE5(end-2:end);
x0 = 2*x4-e4-x5-xd40-xd45;
[ ~, f40 ] = gen0( x0 );
g4=geng(x4);
g5=geng(x5);
g54 = eye(2);
g45 = eye(2);
f54 = [0;0];
f45 = [0;0];
% u54 = f54+g54*u4;
scrLg5=[eye(2) -g54;-g45 2*eye(2)];
scrLg5inv=[2*eye(2) eye(2);eye(2) eye(2)];
scrLg55=[2*eye(2) eye(2)];
scrLg54=[eye(2) eye(2)];
F5 = [f54;f45+f40];
scrF5=g5*scrLg55*F5-g4*scrLg54*F5;
scrG5=g5*scrLg55-g4*scrLg54;
calF5=g5*scrLg55*F5;
calG5=g5*scrLg55;
[~,sig5GRe5,sig5GRx5,sig5GRe4,~]=ADPSEMREMKNN_Basis5(calE5);    

scrLg4=[2*eye(2) -g54;-g45 eye(2)];
scrLg4inv=[eye(2) eye(2);eye(2) 2*eye(2)];
scrLg44=[eye(2) eye(2)];
scrLg45=[eye(2) 2*eye(2)];

scrG4=2*g4*scrLg44-g5*scrLg45;

Gsig5=(scrG5*[eye(2);zeros(2)])'*sig5GRe5'+...
    (scrG4*[zeros(2);eye(2)])'*sig5GRe4'+...
    (calG5*[eye(2);zeros(2)])'*sig5GRx5';
mu5=-0.5*(R5\Gsig5)*Wa5H;
end
