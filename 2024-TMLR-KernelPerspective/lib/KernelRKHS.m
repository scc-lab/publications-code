% The KernelRKHS class implements the following kernels:
%
% 'Exponential with parameter mu: exp(x.'*y/mu)', 
% 'Gaussian with parameter mu: exp(norm(x-y)^2/mu)',
% 'Polynomial with parameters [mu,d,b]: (b + x.'*y/mu)^d', and
% 'Linear with parameter mu = Polynomial with parameters [mu,1,0]'
%
% The kernels functions are written to accept column vector, matrix, 
% or 3D array inputs. The two inputs X and Y have to be the same size in 
% dimensions 1 and 2. If they have different number of pages, one of them 
% has to be a single page. 
% 
% Let the kernel function be denoted by k. The class method obj.K returns
% a u x u x p matrix M, where u = size(X,2), and 
% p = max(size(X,3), size(Y,3)), with 
% 
% M(i,j,k) = k(X(:,i,k),Y(:,j,k)) if if X and Y are same size,
% M(i,j,k) = k(X(:,i,k),Y(:,j,1)) if size(Y,3) = 1, and
% M(i,j,k) = k(X(:,i,1),Y(:,j,k)) if size(X,3) = 1.
%
% The class also includes a method obj.G2K that calculates the inner 
% product of a given vector with the gradient of the kernel function obj.K,
% with respect to the second argument. This function is compatible
% with vector, matrix, or 3D array inputs for the first two arguments. 
% ***The third argument has to be a column vector.***
% 
% Let the gradient $\frac{\partial k(x,y)}{\partial y}$ be denoted by Gk.
% obj.G2K(X,Y,z) returns a matrix N, the same size as M above, with
% 
% N(i,j,k) = z^T * Gk(X(:,i,k),Y(:,j,k)) if X and Y are same size,
% N(i,j,k) = z^T * k(X(:,i,k),Y(:,j,1)) if size(Y,3) = 1, and
% N(i,j,k) = z^T * k(X(:,i,1),Y(:,j,k)) if size(X,3) = 1.
%
% © Rushikesh Kamalapurkar
%
classdef KernelRKHS
    properties
        parameter; % kernel width
        degree; % kernel width
        bias; % kernel width
        type; % kernel type
    end
    methods
        function obj = KernelRKHS(type,parameter)
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
            if isequal(obj.type,'Polynomial')
                if numel(parameter)==1
                    obj.parameter = parameter;
                    obj.degree = 1;
                    obj.bias = 0;
                    warning('Selected linear kernel with parameter mu: k(x,y) = x.''*y/mu')
                elseif numel(parameter)==3
                    obj.parameter = parameter(1);
                    obj.degree = parameter(2);
                    obj.bias = parameter(3);
                else
                    error('The polynomial kernel requires three parameters [mu,d,b], k(x,y) = (b + x.''*y/mu)^d')
                end
            else
                obj.degree = 1;
                obj.bias = 0;
            end
            if isequal(obj.type,'Linear')
                obj.parameter = parameter;
                obj.degree = 1;
                obj.bias = 0;
            end
        end

        function y = K(obj,X,Y)
            % Kernel function
            if isequal(obj.type,'Gaussian')
                y = exp(-1/obj.parameter*(pagetranspose(sum(X.^2,1)) + sum(Y.^2,1) - 2*pagemtimes(X,'transpose',Y,'none')));
            elseif isequal(obj.type,'Exponential')
                y = exp(1/obj.parameter*pagemtimes(X,'transpose',Y,'none'));
            elseif isequal(obj.type,'Linear') || isequal(obj.type,'Polynomial')
                y = (obj.bias + 1/obj.parameter*pagemtimes(X,'transpose',Y,'none')).^obj.degree;
            else
                error(['Kernel type' obj.type 'not implemented']);
            end
        end

        function y = G2K(obj,X,Y,Z)
            % A function that calculates the product of a given vector and
            % the gradient of K with respect to the second argument.
            if isequal(obj.type,'Gaussian')
                y = (2/obj.parameter)*(pagemtimes(X,'transpose',Z,'none') - pagemtimes(Z,'transpose',Y,'none')).*obj.K(X,Y);
            elseif isequal(obj.type,'Exponential')
                y = 1/obj.parameter*pagemtimes(X,'transpose',Z,'none').*obj.K(X,Y);
            elseif isequal(obj.type,'Linear') || isequal(obj.type,'Polynomial')
                if obj.degree == 1
                    y = 1/obj.parameter*pagemtimes(X,'transpose',Z,'none').*ones(size(X,2),size(Y,2),max(size(X,3),size(Y,3)));
                else
                    y = obj.degree/obj.parameter*pagemtimes(X,'transpose',Z,'none').*(obj.bias + 1/obj.parameter*pagemtimes(X,'transpose',Y,'none')).^(obj.degree-1);
                end
            else
                error(['Kernel type' obj.type 'not implemented']);
            end
        end

        function y = OKGram(obj,X,wx,Y,wy)
            y=zeros(size(X,3),size(Y,3));
            for i=1:size(Y,3)
                y(:,i) = squeeze(pagemtimes(pagemtimes(wx(:,1,i).',obj.K(X(:,:,i),Y)),wy));
            end
        end

        function y = OK(obj,x,X,w)
            y=zeros(size(X,3),size(x,2));
            for i=1:size(x,2)
                y(:,i) = squeeze(pagemtimes(obj.K(x,X),w));
            end
        end
    end
end