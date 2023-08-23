% kpcabound(data,sigma,numev,outlier)
%
% kpcabound demonstrates 'Kernel PCA for novelty detection'
%
% kpcabound plots the reconstruction error in feature space into the 
% original space and plots the decision boundary enclosing all data points
% 
% input:
%
% data: array of data points, one row for each data point
% sigma: width of Gaussian kernel
% numev: number of eigenvalues to be extracted
% outlier: number of points outside the decision boundary
%
% (c) Heiko Hoffmann
%
% This code is distributed in the hope that it will be useful, but 
% without any warranty; even without the implied warranty of 
% merchantability or fitness for any purpose.
%
% Slight modification by Rushikesh Kamalapurkar

function [RTrain,RTest] = kpcabound(data,sigma,numev,test)

     
 [n,d] = size(data);
 % n : number of data points
 % d : dimension of data points
 
 % kernel matrix:
 K = zeros(n,n);
 
 % kernel parameter:
 param = 0.5/(sigma*sigma);
 
 %fprintf('computing kernel matrix K\n');
 for i=1:n 
   for j=i:n
     K(i,j) = kernel(data(i,:),data(j,:),param);
     K(j,i) = K(i,j);
   end
 end
 
 % correct K for non-zero center of data in feature space:
 Krow = sum(K,1)/n;
 Ksum = sum(Krow)/n;
 
 for i=1:n
   for j=1:n
     K(i,j) = K(i,j) - Krow(i) - Krow(j) + Ksum;
   end
 end
 
 %fprintf('extracting eigenvectors of K\n');
 opts.disp = 0;
 [alpha,lambda]=eigs(K,numev,'lm',opts);
 
 % residual variance:
 resvar = (trace(K)-trace(lambda));
 %fprintf('residual variance relative to total variance in feature space: %f %%\n',100*resvar/trace(K));
	 
 % normalize alpha:
 alpha = alpha * inv(sqrt(lambda));

 % compute some helper vectors:
 sumalpha = sum(alpha,1);
 alphaKrow = Krow * alpha;
 
 %fprintf('evaluating reconstruction error for all data points\n');
 
 RTrain = zeros(n,1);
 for i=1:n
   x = data(i,:); % test point
   RTrain(i) = recerr(x,data,param,alpha,alphaKrow,sumalpha,Ksum);
 end
 [m,~] = size(test); 
 RTest = zeros(m,1);
 for i=1:m
   x = test(i,:); % test point
   RTest(i) = recerr(x,data,param,alpha,alphaKrow,sumalpha,Ksum);
 end
end