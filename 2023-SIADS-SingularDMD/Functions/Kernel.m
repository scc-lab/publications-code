% Code written by Rushikesh Kamalapurkar
% The Kernel class implements the following kernels:
%
% 'Exponential', 
% 'Gaussian'
% 'Linear'
%
% The kernels functions are written to accept column vector, matrix, 
% or 3D array inputs. If the kernel function is k, then the class function
% Kernel.K returns the following outputs based on the size of the inputs
% 
% n x 1 column vector inputs:
%   Kernel.K(x,y) returns the scalar k(x,y)
% n x m matrix inputs X = [x1 x2 ... xm], Y = [y1 y2 ... ym]:
%   Kernel.K(X,Y) treats the columns of the input matrices as vectors and 
%   returns the m x m matrix M(i,j) = k(xi,yj)
% n x m x p 3D array inputs X and Y:
%   Kernel.K(X,Y) operates on each page of X and Y as above to generate the
%   m x m x p array M(i,j,k) = k(xi_k,yj_k)
%   where xi_k is the i-th column on the k-th page of X.
classdef Kernel
    properties
        parameter; % kernel width
        type; % kernel type
        K; % kernel function
    end
    methods
        function obj = Kernel(type,parameter)
            % Processing arguments and setting defaults
            if nargin == 2
                obj.type = type;
                obj.parameter = parameter;
            elseif nargin == 0
                obj.type = 'Gaussian'; % default
                obj.parameter = 1; % default
                warning('Setting kernel type to Gaussian and kernel parameter to 1');
            else
                error('Two (or zero) input arguments needed');
            end
            % Kernel functions
            if isequal(type,'Gaussian')
                obj.K = @(X,Y) exp(-1/parameter*(pagetranspose(sum(X.^2,1)) + sum(Y.^2,1) - 2*pagemtimes(X,'transpose',Y,'none')));
            elseif isequal(type,'Exponential')
                obj.K = @(X,Y) exp(1/parameter*pagemtimes(X,'transpose',Y,'none'));
            elseif isequal(type,'Linear')
                obj.K = @(X,Y) 1/parameter*pagemtimes(X,'transpose',Y,'none');
            else
                error(['Kernel type' type 'not implemented']);
            end
        end
    end
end
                
