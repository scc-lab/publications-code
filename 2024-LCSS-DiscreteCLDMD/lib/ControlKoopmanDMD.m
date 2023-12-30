% This function performs dynamic mode decomposition of a controlled 
% dynamical system from snapshots (x,u,y) of the form
% 
% y = f(x) + g(x)u
%
% [Z,L,ef,cr,dr,fc] = controlKoopmanDMD(KT,K,X,U,Y,deltaT,mu) OR
% [Z,L,ef,cr,dr,fc] = controlKoopmanDMD(KT,K,X,U,Y,deltaT,mu,l)
%
% Inputs:
%    1) KT:      svRKHS kernel 
%    2) K:       vvRKHS kernel 
%    3) X:       Input data
%                First dimension: State
%                Second dimension: Snapshot
%    4) U:       Control input data
%                First dimension: Control
%                Second dimension: Snapshot
%    5) Y:       Output data, Y(:,k) = f(X(:,k)) + g(x)U(:,k)
%                First dimension: State
%                Second dimension: Snapshot
%    6) deltaT:  Sample time
%    6) mu:      A feedback control function compatible with matrix inputs
%    7) l:       (Optional argument, default=0) Regularization coefficient
%                (needed if the Gram matrix is rank deficient)
%
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
function [Z,L,ef,cr,dr,fc] = controlKoopmanDMD(KT,K,X,U,Y,deltaT,mu,varargin)

% Processing optional arguments and setting defaults
if nargin == 7
    l = 0; % default
elseif nargin == 8
    l = varargin{1};
elseif nargin > 8
    error("Too many input arguments.")
end

% Values of the feedback controller
MU = mu(X);

U = permute([ones(1,size(U,2));U],[2 1 3]);
MU = permute([ones(1,size(MU,2));MU],[2 1 3]);

% vvRKHS Gram matrix and interaction matrix - M x M
G = sum(squeeze(K.K(X,X)).*pagemtimes(U,'none',U,'transpose'),3);
I = sum(squeeze(K.K(X,X)).*pagemtimes(MU,'none',U,'transpose'),3);

% svRKHS Gram matrix and interaction matrix - M x M
GT = KT.K(X,X);
IT = KT.K(X,Y);

% DMD
GT=GT+l*eye(size(GT)); % Regularization
G=G+l*eye(size(G)); % Regularization
[V,L] = eig(GT\(I*(G\IT.'))); % Eigendecomposition of the finite rank representation
contL = log(diag(L))./deltaT; % Continuous-time eigenvalues
C = V./diag(sqrt(V'*GT*V)).'; % Normalized eigenvectors of finite rank representation
Z = X/(C.'*GT); % Koopman modes

% Reconstruction
ef = @(x) (KT.K(x,X)*C).';% Koopman eigenfunctions
% Discrete Koopman reconstruction function (or vector field if k = 1)
dr = @(k,x) real(Z*(ef(x).*diag(L).^k));
% Continuous Koopman reconstruction function
cr = @(t,x) real(Z*(ef(x).*exp(contL*t)));
% Continuous Koopman vector field
fc = @(x) real(Z*(ef(x).*contL));
end