% This is an implementation of the Koopman DMD method in 
% 
% https://arxiv.org/abs/2106.00106
%
% Code written by Rushikesh Kamalapurkar

function [Z,L,ef,cr,dr] = KoopmanDMD(X,Y,K,deltaT,varargin)
% This function computes the modes, eigenvalues, and eigenfunction
% coefficients of the Koopman operator that relates the inputs X to the
% outputs Y
% [Z,L,ef,cr,dr] = WilliamsKDMD(X,Y,K,deltaT) OR
% [Z,L,ef,cr,dr] = WilliamsKDMD(X,Y,K,deltaT,l)
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
%          Example 3, linear dot product
%             K.K = @(X,Y) 1/mu*pagemtimes(X,'transpose',Y,'none');
%       deltaT: Sample time
%       l: (optional argument, default = 0) Regularization coefficient
%          (needed if the Gram matrix is rank deficient)
% Outputs:
%       Z: Koopman modes
%       L: Koopman eigenvalues
%       ef: Koopman eigenfunctions
%       cr: Continuous Koopman reconstruction function, cr(t,x0)
%           reconstructs the system state at time t starting from x0
%       dr: Discrete Koopman reconstruction function, dr(k,x0) reconstructs
%           the k-th snapshot starting from x0

% Processing optional arguments and setting defaults 
if nargin == 4
    l = 0; % default
elseif nargin == 5
    l = varargin{1};
elseif nargin > 5
    error("Too many input arguments.")
end
% Gram matrix and interaction matrix
G = K.K(X,X);
I = K.K(X,Y);

% DMD
G = G + l*eye(size(G));
[V,L] = eig(G\I'); % Eigendecomposition of finite rank representation
contL = log(diag(L))./deltaT;
C = V./diag(sqrt(V'*G*V)).'; % Koopman eigenfunction coefficients
Z = X/(C.'*G); % Koopman modes

% Koopman eigenfunctions
ef = @(x0) (K.K(x0,X)*C).';

% Koopman reconstruction functions
dr = @(k,x0) Z*(ef(x0).*diag(L).^k);
cr = @(t,x0) Z*(ef(x0).*exp(contL*t));
end
