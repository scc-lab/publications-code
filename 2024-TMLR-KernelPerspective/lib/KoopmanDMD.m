% This function computes the modes, eigenvalues, and eigenfunction
% coefficients of a finite-rank representation of the Koopman operator
% for the vector field that relates the inputs X to the outputs Y. This is 
% an efficient implementation of the kernel DMD method in
%
% https://arxiv.org/abs/2106.00106
%
% [Z,L,ef,cr,dr,fc] = KoopmanDMD(X,Y,K,deltaT) OR
% [Z,L,ef,cr,dr,fc] = KoopmanDMD(X,Y,K,deltaT,l)
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
%          Example 3, linear dot product
%             K.K = @(X,Y) 1/mu*pagemtimes(X,'transpose',Y,'none');
%       deltaT: Sample time
%       l: (optional argument, default = 0) Regularization coefficient
%          (needed if the Gram matrix is rank deficient)
% Outputs:
%       Z:  Koopman modes
%       L:  Koopman eigenvalues
%       ef: Koopman eigenfunctions
%       cr: Continuous Koopman reconstruction function, cr(t,x0)
%           reconstructs the system state at time t starting from x0
%       dr: Discrete Koopman reconstruction function, dr(k,x0) reconstructs
%           the k-th snapshot starting from x0. With k = 1, we get the
%           vector field $ f_{d} $ of the estimated discrete time system as
%           $ x_{k+1} = fd(x_{k}) = dr(1,x_{k}) $
%       fc: The vector field of the estimated continuous time system 
%           $ \dot{x} = fc(x) $. This is relevant only if the data were
%           sampled from a continuous time system.
%
% © Rushikesh Kamalapurkar
%
function [Z,L,ef,cr,dr,fc] = KoopmanDMD(X,Y,K,deltaT,varargin)
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
    % DMD (see https://arxiv.org/abs/2106.00106 for details)
    G = G + l*eye(size(G)); % Regularization of the Gram matrix
    [V,L] = eig(G\I'); % Eigendecomposition of finite rank representation
    contL = log(diag(L))./deltaT; % Continuous-time eigenvalues
    C = V./diag(sqrt(V'*G*V)).'; % Koopman eigenfunction coefficients
    Z = X/(C.'*G); % Koopman modes
    ef = @(x0) (K.K(x0,X)*C).';% Koopman eigenfunctions
    % Discrete Koopman reconstruction function (or vector field if k = 1)
    dr = @(k,x0) real(Z*(ef(x0).*diag(L).^k));
    % Continuous Koopman reconstruction function
    cr = @(t,x0) real(Z*(ef(x0).*exp(contL*t)));
    % Continuous Koopman vector field
    fc = @(x0) real(Z*(ef(x0).*contL));
end
