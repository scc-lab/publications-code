% This is an implementation of the kernel DMD method in 
% Williams et al., 2015, see
% https://doi.org/10.3934/jcd.2015005
%
% Code written by Rushikesh Kamalapurkar

function [Z,L,ef,cr,dr] = WilliamsKDMD(X,Y,K,deltaT,varargin)
% This function computes the modes, eigenvalues, and eigenfunction
% coefficients of a finite-rank representation of the Koopman operator
% for the vector field that relates the inputs X to the outputs Y
%
% [Z,L,ef,cr,dr] = WilliamsKDMD(X,Y,K,deltaT)
%
% Inputs:
%       X: Input data,
%           First dimension: state variables
%           Second dimension: snapshots
%       Y: Output data, Y = f(X)
%       K: Object of the class 'Kernel', where K.K is the kernel function
%          Example 1, Gaussian:
%             K.K = @(X,Y) exp(-1/mu*(pagetranspose(sum(X.^2,1)) + ...
%             sum(Y.^2,1) - 2*pagemtimes(X,'transpose',Y,'none'))); 
%          Example 2, exponential dot product:
%             K.K = @(X,Y) exp(1/mu*pagemtimes(X,'transpose',Y,'none'));
%       deltaT: Sample time
% Outputs:
%       Z: Koopman modes
%       L: Koopman eigenvalues
%       ef: Koopman eigenfunctions
%       cr: Continuous Koopman reconstruction function, cr(t,x0)
%           reconstructs the system state at time t starting from x0
%       dr: Discrete Koopman reconstruction function, dr(k,x0) reconstructs
%           the k-th snapshot starting from x0

% Gram matrix and interaction matrix
GHat = K.K(X,X);
AHat = K.K(X,Y).';

[Q,SigmaSquared] = eig(GHat);
Sigma = sqrt(SigmaSquared);
KHat = (pinv(Sigma)*Q')*AHat*(Q*pinv(Sigma));
[VHat,L] = eig(KHat);
mu = log(diag(L))./deltaT;
ef = @(x) K.K(x,X)*(Q*pinv(Sigma)*VHat); % eigenfunction at x
Z = (pinv(VHat)*pinv(Sigma)*Q'*X').'; % Koopman modes

% Koopman reconstruction function
cr = @(t,x0) Z*(ef(x0).'.*exp(mu*t));
dr = @(k,x0) Z*(ef(x0).'.*diag(L).^k);
end
