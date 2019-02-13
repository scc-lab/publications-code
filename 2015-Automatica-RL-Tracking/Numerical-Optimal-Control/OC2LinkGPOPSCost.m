function [M,L] = OC2LinkGPOPSCost(solcost)

Q=[10  0   0   0; 
   0   10  0   0;
   0   0   2   0;
   0   0   0   2]; % State Penalty
R = [1 0;0 1]; % Control Penalty
zeta = solcost.state;
e = zeta(:,1:4);
mu = solcost.control;
% Calculate cost. 
L = (dot(e,e*Q.',2)+dot(mu,mu*R.',2)); % Note transposes
M = 0;

