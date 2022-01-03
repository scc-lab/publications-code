function [Z,L,C,f] = LiouvilleDMD(K,W,t,varargin)
% This function performs dynamic mode decomposition of a dynamical system
% from sampled trajectories of the system. Trajectories can be of different
% lengths and irregularly sampled, but
% *the number of samples in each trajectory needs to be odd.*
%
% [Z,L,C,f] = LiouvilleDMD(K,W,T) OR
% [Z,L,C,f] = LiouvilleDMD(K,W,T,a) OR
% [Z,L,C,f] = LiouvilleDMD(K,W,T,a,l)
%
% Inputs:
%    1) K: A kernel function (compatible with 3D arrays)
%          Example 1, Gaussian:
%             K = @(X,Y) exp(-1/mu*(pagetranspose(sum(X.^2,1)) + ...
%             sum(Y.^2,1) - 2*pagemtimes(X,'transpose',Y,'none'))); 
%          Example 2, exponential dot product:
%             K = @(X,Y) exp(1/mu*pagemtimes(X,'transpose',Y,'none'));
%    2) W: A dataset of trajectories (3D array)
%          First dimension: State 
%          Second dimension: Time (size = length of longest trajectory)
%          Third dimension: Trajectory number
%    3) t: A matrix of sample times
%          First dimension: Time (size = length of longest trajectory)
%          Second dimension: Trajectory number
%    4) a: (Optional argument, default=1) Scaling factor 
%    5) l: (Optional argument, default=0) Regularization coefficient
%
% *Shorter trajectories need to be padded with zeros.*
% *Sample times for shorter trajectories need to be padded with NaNs.*
%
% Outputs:
%    1) Z: Liouville modes (State dimension x number of modes)
%    2) L: Eigenvalues (number of modes x 1)
%    3) C: Eigenfunction coefficients (number of modes x number of modes)
%    4) f: Reconstruction function

if nargin == 3
    a = 1;
    l = 0;
elseif nargin == 4
    a = varargin{1};
    l = 0;
elseif nargin == 5
    a = varargin{1};
    l = varargin{2};
elseif nargin > 5
    error("Too many input arguments.")
end
    
N = size(W,3); % Total number of trajectories

% Store trajectory lengths for interaction matrix calculation
Lengths = size(t,1)-sum(isnan(t));

% Simpsons rule weights
S = reshape(genSimpsonsRuleWeights(t,1),size(t,1),1,size(t,2));

% Gram matrix and interaction matrix
G=zeros(N);
I=zeros(N);
for i=1:N
    G(:,i) = squeeze(pagemtimes(pagemtimes(S(:,1,i).',K(W(:,:,i),W)),S));
    I(:,i) = pagemtimes(K(a*W(:,Lengths(i),i),W) - K(a*W(:,1,i),W),S);
end

% DMD
G = G + l*eye(size(G)); % Regularization
[V,D] = eig(G\I.'); % Eigendecomposition of finite rank representation
C = V./diag(sqrt(V'*G*V)).'; % Liouville eigenfunction coefficients
IntMat = squeeze(pagemtimes(W,S)); % Integrals of trajectories
Z = IntMat/(C.'*G); % Liouville modes
L = diag(D); % Liouville eigenvalues

% Reconstruction
% Occupation kernels evaluated at x0: squeeze(pagemtimes(K(x0,W),S))
% Eigenfunctions evaluated at x0: C.'*squeeze(pagemtimes(K(x0,W),S))
% Reconstruction function:
f = @(t,x0) Z*((C.'*squeeze(pagemtimes(K(x0,W),S))).*exp(L*t)); 

end