% This function performs dynamic mode decomposition of a second order 
% dynamical system from sampled trajectories. For details, see
% 
% https://arxiv.org/abs/2101.02646
%
% [Z,D,ef,r,f] = SecondOrderLiouvilleDMD(K,G2K,X,t) OR
% [Z,D,ef,r,f] = SecondOrderLiouvilleDMD(K,G2K,X,t,XDot0) OR
% [Z,D,ef,r,f] = SecondOrderLiouvilleDMD(K,G2K,X,t,XDot0,l)
%
% Inputs:
%    1) K: An object of the class Kernel.
%          K.K needs to be the kernel function, compatible with 3D arrays. 
%          K.G2K needs to be a function that calculates the product of a 
%          given vector and the gradient of K with respect to the second
%          argument.
%          Example 1, Gaussian:
%             K.K = @(X,Y) exp(-1/mu*(pagetranspose(sum(X.^2,1)) + ...
%             sum(Y.^2,1) - 2*pagemtimes(X,'transpose',Y,'none'))); 
%             K.G2K = @(X,Y,Z) pagemtimes(2*(X - Y)/mu,...
%              'transpose', Z, 'none').*K(X,Y);
%          Example 2, exponential dot product:
%             K.K = @(X,Y) exp(1/mu*pagemtimes(X,'transpose',Y,'none'));
%             K.G2K = @(X,Y,Z) 1/mu * ...
%              pagemtimes(X,'transpose',Z,'none').*K(X,Y);
%          Example 2, linear dot product:
%             K.K = @(X,Y) 1/mu*pagemtimes(X,'transpose',Y,'none');
%             K.G2K = @(X,Y,Z) 1/mu*pagemtimes(X,'transpose',Z,'none');
%    2) X: A dataset of trajectories (3D array)
%          First dimension: State 
%          Second dimension: Time (size = length of longest trajectory)
%          Third dimension: Trajectory number
%
% *Trajectories can be of different lengths and irregularly sampled.*
% *the number of samples in each trajectory needs to be odd.*
% *Shorter trajectories need to be padded with zeros.*
%
%    3) t: A matrix of sample times
%          First dimension: Time (size = length of longest trajectory)
%          Second dimension: Trajectory number
%
% *Sample times for shorter trajectories need to be padded with NaNs.*
%
%    5) XDot0: (Optional argument) Initial time derivatives
%          First dimension: State 
%          Second dimension: Trajectory number
%    6) l: (Optional argument) Regularization coefficient
%
% Outputs:
%          1) Z: Liouville modes
%          2) D: Diagonal matrix of Liouville eigenvalues
%          3) ef: Liouville eigenfunctions
%          4) r: Trajectory prediction function
%          5) f: Vector Field
%
% © Rushikesh Kamalapurkar, Ben Russo, and Joel Rosenfeld
%
function [Z,D,ef,r,f] = SecondOrderLiouvilleDMD(K,X,t,varargin)
if nargin == 3
    if size(X,1)==1
        XDiff = (squeeze(X(:,2,:) - X(:,1,:))).';
    else
        XDiff = squeeze(X(:,2,:) - X(:,1,:));
    end
    XDot0 = XDiff./(t(2,:) - t(1,:));
    l = 0;
elseif nargin == 4
    XDot0 = varargin{1};
    l = 0;
elseif nargin == 5
    XDot0 = varargin{1};
    if isempty(XDot0)
        if size(X,1)==1
            XDiff = (squeeze(X(:,2,:) - X(:,1,:))).';
        else
            XDiff = squeeze(X(:,2,:) - X(:,1,:));
        end
        XDot0 = XDiff./(t(2,:) - t(1,:));
    end
    l = varargin{2};
elseif nargin > 5
    error("Too many input arguments");
end
    
M = size(X,3); % Total number of trajectories
n = size(X,1); % State Dimension
N = size(t,1) - sum(isnan(t)); % Trajectory lengths

% Simpsons rule weights
w = reshape(genSimpsonsRuleWeights(t,1),size(t,1),1,size(t,2));
t(isnan(t))=0;
T=zeros(1,1,size(t,2));
for i=1:size(t,2)
    if N(i) == 1
        T(1,1,i) = t(1,i);
    else
        T(1,1,i) = t(find(t(:,i),1,'last'),i);
    end
end
t = reshape(t,size(t,1),1,size(t,2));

% Gram matrix and interaction matrix
G=zeros(M);
I=zeros(M);
wT=w.*(T-t);
for j=1:M
    G(:,j) = squeeze(pagemtimes(wT,'transpose',pagemtimes(K.K(X,X(:,:,j)),wT(:,1,j)),'none'));
    I(:,j) = squeeze(pagemtimes(wT,'transpose',K.K(X,X(:,N(j),j)) - K.K(X,X(:,1,j)) - T(1,1,j).*K.G2K(X,X(:,1,j),XDot0(:,j)),'none'));
end

% DMD
G = G + l*eye(size(G)); % Regularization
[V,D] = eig(G\I.'); % Eigendecomposition of finite rank representation
C = V./diag(sqrt(V'*G*V)).'; % Liouville eigenfunction coefficients
IntMat = reshape(pagemtimes(X,wT),n,M); % Integrals of trajectories
Z = IntMat/(C.'*G); % Liouville modes
ef = @(x) C.'*squeeze(pagemtimes(wT,'transpose',K.K(X,x),'none')); % Eigenfunctions evaluated at x
% efDot = @(x) C.'*squeeze(pagemtimes(wT,'transpose',K.G2K(X,x0,xDot0),'none')); % Time derivative of eigenfunctions at t = 0

% SysID

r = @(t,x0,xDot0) real(Z*(... % Trajectory prediction
    0.5*(C.'*squeeze(pagemtimes(wT,'transpose',K.K(X,x0),'none'))...
    + (C.'*squeeze(pagemtimes(wT,'transpose',K.G2K(X,x0,xDot0),'none')))./sqrt(diag(D))).*exp(sqrt(diag(D))*t) +...
    0.5*(C.'*squeeze(pagemtimes(wT,'transpose',K.K(X,x0),'none'))...
    - (C.'*squeeze(pagemtimes(wT,'transpose',K.G2K(X,x0,xDot0),'none')))./sqrt(diag(D))).*exp(-sqrt(diag(D))*t)));
f = @(x) real(Z*D*C.'*squeeze(pagemtimes(wT,'transpose',K.K(X,x),'none'))); % Vector field
end
