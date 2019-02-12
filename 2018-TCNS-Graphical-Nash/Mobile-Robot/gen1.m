function [ mu1,scrG1,scrF1,calG1,calF1,Gsig1 ] = gen1( calE1,Wa1H,R1 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
e1 = calE1(1:3);
x1 = calE1(4:end);
x0 = x1-e1;
[ f0, f10 ] = gen0( x0 );
g1=[cos(x1(3)) 0;sin(x1(3)) 0;0 1];
F1 = f10;
scrF1=g1*F1-f0;
scrG1=g1;
calF1=g1*F1;
calG1=g1; 
[~,sig1GRe1,sig1GRx1,~]=ADPSEMREMKNN_Basis1(calE1);         
Gsig1=scrG1'*sig1GRe1'+calG1'*sig1GRx1';
mu1=-0.5*(R1\Gsig1)*Wa1H;
end

