% This is an implementation of occupation kernel PCA for fault detection'
%
% [Rte,RTr] = OKPCAReconstructionError(K,Tr,tTr,Te,tTe,N)
%
% Inputs:
%   (1) K: Object of the class 'Kernel', where K.K is the kernel function
%          Example 1, Gaussian:
%             K.K = @(X,Y) exp(-1/mu*(pagetranspose(sum(X.^2,1)) + ...
%             sum(Y.^2,1) - 2*pagemtimes(X,'transpose',Y,'none'))); 
%          Example 2, exponential dot product:
%             K.K = @(X,Y) exp(1/mu*pagemtimes(X,'transpose',Y,'none'));
%
% (2,4) Tr,Te: Datasets of training and test trajectories (3D array)
%          First dimension: State 
%          Second dimension: Time (size = length of longest trajectory)
%          Third dimension: Trajectory number
%
% *Trajectories can be of different lengths and irregularly sampled.*
% *The number of samples in each trajectory needs to be odd.*
% *Shorter trajectories need to be padded with zeros.*
%
% (3,5) tTr,tTe: Matrices of sample times for training and test trajectories
%          First dimension: Time (size = length of longest trajectory)
%          Second dimension: Trajectory number
%
% *Sample times for shorter trajectories need to be padded with NaNs.*
%
%   (6) N: Number of principle components to use
%
% Outputs:
%    1) RTest: Reconstruction error for test data
%    2) RTrain: Reconstruction error for training data
%
% © Rushikesh Kamalapurkar
%
function [Rte,RTr] = OKPCAReconstructionError(K,Tr,tTr,Te,tTe,N)
    M = size(Tr,3);
    P = size(Te,3);

    % Simpsons rule weights
    STr = reshape(genSimpsonsRuleWeights(tTr,1),size(tTr,1),1,size(tTr,2));
    STe = reshape(genSimpsonsRuleWeights(tTe,1),size(tTe,1),1,size(tTe,2));
    
    % Gram matrix for training data 
    G=zeros(M);
    % Inner product matrix for training and testing data
    teTr = zeros(P,M);
    for j=1:M
        G(:,j) = squeeze(pagemtimes(pagemtimes(STr(:,1,j).',K.K(Tr(:,:,j),Tr)),STr));
        teTr(:,j) = pagemtimes(pagemtimes(STe,'transpose',K.K(Te,Tr(:,:,j)),'none'),STr(:,:,j));
    end
    
    % Diagonal of the Gram matrix for testing data
    teTe = zeros(P,1);
    for i=1:P
        teTe(i,1) = STe(:,1,i).'*K.K(Te(:,:,i),Te(:,:,i))*STe(:,1,i);
    end

    % Occupation kernel principle component analysis
    % Data centering
    JM = (1/M)*ones(M);
    KCentered = G - JM*G - G*JM + JM*G*JM;

    % Eigendecomposition
    [alpha,lambda]=eigs(KCentered,N,'largestabs');

    % Normalization of eigenvectors
    alpha = alpha/sqrt(lambda);
    
    % Reconstruction error for testing data
    normPhiTildeSq = teTe - (2/M)*sum(teTr,2) + (1/M^2)*sum(G,"all");
    PhiTildev = (teTr - (1/M)*sum(G,1) - (1/M)*sum(teTr,2) + (1/M^2)*sum(G,"all"))*alpha;
    Rte = normPhiTildeSq - sum(PhiTildev.^2,2);
    
    % Reconstruction error for training data
    normPhiTildeSqTr = diag(G) - (2/M)*sum(G,2) + (1/M^2)*sum(G,"all");
    PhiTildevTr = (G - (1/M)*sum(G,1) - (1/M)*sum(G,2) + (1/M^2)*sum(G,"all"))*alpha;
    RTr = normPhiTildeSqTr - sum(PhiTildevTr.^2,2);
end