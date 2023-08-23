% k=kernel(x,y,param)
%
% Kernel function for KPCA. Here, a Gaussian kernel is implemented.
% You may changed this function to try other kernels.
%
% x,y: two data points
% param: parameter of kernel function

function k=kernel(x,y,param)
  diff = x-y;
  k = exp(-(diff * diff')*param);
