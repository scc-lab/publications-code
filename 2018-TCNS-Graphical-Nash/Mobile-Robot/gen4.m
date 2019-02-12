function [ mu4,scrG4,scrF4,calG4,calF4,Gsig4,F4,x0 ] = gen4( calE4,calE5,Wa4H,R4,xd40,xd45 )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
e4 = calE4(1:3);
x4 = calE4(end-2:end);
x5 = calE5(end-2:end);
x0 = 2*x4-e4-x5-xd40-xd45;
[ f0, f40 ] = gen0( x0 );
g5=geng(x5);
g4=geng(x4);
g45 = eye(2);
g54 = eye(2);
f45 = [0;0];
f54 = [0;0];
% u45 = f45+g45*u5;

scrLg5=[eye(2) -g45;-g54 2*eye(2)];
scrLg5inv=[2*eye(2) eye(2);eye(2) eye(2)];
scrLg4=[2*eye(2) -g45;-g54 eye(2)];
scrLg4inv=[eye(2) eye(2);eye(2) 2*eye(2)];
scrLg55=[2*eye(2) eye(2)];
scrLg54=[eye(2) eye(2)];
scrLg44=[eye(2) eye(2)];
scrLg45=[eye(2) 2*eye(2)];

F4 = [f40+f45;f54];

scrF4=2*g4*scrLg44*F4-f0-g5*scrLg45*F4;
scrG4=2*g4*scrLg44-g5*scrLg45;
calF4=g4*scrLg44*F4;
calG4=g4*scrLg44;
[~,sig4GRe4,sig4GRx4,sig4GRe5,~]=ADPSEMREMKNN_Basis4(calE4);    

scrG5=g5*scrLg55-g4*scrLg54;

Gsig4=(scrG4*[eye(2);zeros(2)])'*sig4GRe4'+...
    (scrG5*[zeros(2);eye(2)])'*sig4GRe5'+...
    (calG4*[eye(2);zeros(2)])'*sig4GRx4';
mu4=-0.4*(R4\Gsig4)*Wa4H;