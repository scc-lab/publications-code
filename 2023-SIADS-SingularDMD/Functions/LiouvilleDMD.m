% This is an implementation of the Liouville DMD method in 
% Rosenfeld et al. 2022, see
% https://link.springer.com/article/10.1007/s00332-021-09746-w
%
% Code written by Rushikesh Kamalapurkar

function [Z,L,ef,r,f] = LiouvilleDMD(K,W,t,varargin)
% This function performs dynamic mode decomposition of a dynamical system
% from sampled trajectories of the system. 
%
% [Z,L,ef,r,f] = LiouvilleDMD(K,W,T) OR
% [Z,L,ef,r,f] = LiouvilleDMD(K,W,T,a) OR
% [Z,L,ef,r,f] = LiouvilleDMD(K,W,T,a,l)
%
% Inputs:
%    1) K: Object of the class 'Kernel', where K.K is the kernel function
%          Example 1, Gaussian:
%             K.K = @(X,Y) exp(-1/mu*(pagetranspose(sum(X.^2,1)) + ...
%             sum(Y.^2,1) - 2*pagemtimes(X,'transpose',Y,'none'))); 
%          Example 2, exponential dot product:
%             K.K = @(X,Y) exp(1/mu*pagemtimes(X,'transpose',Y,'none'));
%
%    2) W: A dataset of trajectories (3D array)
%          First dimension: State 
%          Second dimension: Time (size = length of longest trajectory)
%          Third dimension: Trajectory number
%
% *Trajectories can be of different lengths and irregularly sampled.*
% *The number of samples in each trajectory needs to be odd.*
% *Shorter trajectories need to be padded with zeros.*
%
%    3) t: A matrix of sample times
%          First dimension: Time (size = length of longest trajectory)
%          Second dimension: Trajectory number
%
% *Sample times for shorter trajectories need to be padded with NaNs.*
%
%    4) a: (Optional argument, default=1) Scaling factor if using scaled
%          Liouville operators
%
%    5) l: (Optional argument, default=0) Regularization coefficient
%          (needed if the Gram matrix is rank deficient)
%
% Outputs:
%    1) Z: Liouville modes (State dimension x number of modes)
%    2) L: Eigenvalues (number of modes x 1)
%    3) ef: Eigenfunctions
%    4) r: Reconstruction function
%    5) f: Vector field

% Processing optional arguments and setting defaults
if nargin == 3
    a = 1; % default
    l = 0; % default
elseif nargin == 4
    if isempty(varargin{1})
        a = 1;
    else
        a = varargin{1};
    end
    l = 0;
elseif nargin == 5
    if isempty(varargin{1})
        a = 1;
    else
        a = varargin{1};
    end
    l = varargin{2};
elseif nargin > 5
    error("Too many input arguments.")
end
    
N = size(W,3); % Total number of trajectories
n = size(W,1); % State Dimension

% Store trajectory lengths for interaction matrix calculation
Lengths = size(t,1)-sum(isnan(t));

% Simpsons rule weights
S = reshape(genSimpsonsRuleWeights(t,1),size(t,1),1,size(t,2));

% Gram matrix and interaction matrix
G=zeros(N);
I=zeros(N);
for i=1:N
    G(:,i) = squeeze(pagemtimes(pagemtimes(S(:,1,i).',K.K(W(:,:,i),W)),S));
    I(:,i) = pagemtimes(K.K(a*W(:,Lengths(i),i),W) - K.K(a*W(:,1,i),W),S);
end

% DMD
G = G + l*eye(size(G)); % Regularization
[V,D] = eig(G\I.'); % Eigendecomposition of finite rank representation
C = V./diag(sqrt(V'*G*V)).'; % Normalized eigenvectors of finite rank representation
IntMat = reshape(pagemtimes(W,S),n,N); % Integrals of trajectories
Z = IntMat/(C.'*G); % Liouville modes
L = diag(D); % Eigenvalues of the finite-rank representation

% Reconstruction
% Occupation kernels evaluated at x0: squeeze(pagemtimes(K(x0,W),S))
% Eigenfunctions evaluated at x
ef = @(x) C.'*squeeze(pagemtimes(K.K(x,W),S));
% Reconstruction function:
r = @(t,x0) Z*((C.'*squeeze(pagemtimes(K.K(x0,W),S))).*exp(L*t)); 
% Vector field:
f = @(x) real(Z*((C.'*squeeze(pagemtimes(K.K(x,W),S))).*L));
end