% err = recerr(x,data,param,alpha,alphaKrow,sumalpha,Ksum)
% 
% This function computes the reconstruction error of x in feature
% space.
%
% x: test point
% data: data set from which KPCA was computed
% param: kernel parameter
% alpha,alphaKrow,sumalpha,Ksum: resulting data from KPCA
% (see kpcabound.m)

function err = recerr(x,data,param,alpha,alphaKrow,sumalpha,Ksum)
  
   n = size(data,1);
   k = zeros(1,n);
   for j=1:n
     k(j) = kernel(x,data(j,:),param);
   end
   
   % projections:
   f = k*alpha - sumalpha * (sum(k)/n - Ksum) - alphaKrow;
   
   % reconstruction error:
   err = kernel(x,x,param) - 2*sum(k)/n + Ksum - f*f';
