function [ f0, f10 ] = gen0( x0 )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
TimePeriod = 6*pi;
FreqParam = 0.1;
xdD3 = (FreqParam/TimePeriod)*(TimePeriod^2-x0(3)^2);
f0 = [xdD3*cos(x0(3)); xdD3*sin(x0(3)); xdD3];
g0plus=[cos(x0(3)) 0;sin(x0(3)) 0;0 1].'; % the third component of xdij is zero for all i,j
f10 = g0plus*(f0);
end

