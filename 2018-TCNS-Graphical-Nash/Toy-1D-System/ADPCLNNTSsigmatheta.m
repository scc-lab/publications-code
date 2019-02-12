function [sigthi,p_thetai] = ADPCLNNTSsigmatheta(i,xi,t)
if i==0
%     f0=[xi(2)-xi(1); xi(2)-2*xi(1)];
	f0=-0.1*xi;
    %f0=4*sqrt(1-(xi/2)^2);
	%f0=4*sqrt(1-(xi/2)^2);
    sigthi=f0;
else
%     sigthi=[xi(1);xi(2);xi(2)*(1-(cos(2*xi(1))+2)^2)];
    sigthi=[xi^3;xi^2];
end

p_thetai = length(sigthi);