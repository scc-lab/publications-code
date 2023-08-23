% Generate weights for Simpson's rule integration for a function f(t) 
% sampled unevenly at sampling points t(1), t(2), ...
% Need odd number of samples
% Inputs:
%    genSimpsonsRuleWeights(t): A vector containing sampling locations
%    genSimpsonsRuleWeights(t,dim):
%        t: A matrix containing sampling locations, padded with NaN for
%           variable length trajectories
%        dim: dimension of the time increments
% Output: Weights corresponding to the sampling locations, the same size as
% the input
% The code is based on the Composite Simpson's rule for irregularly spaced
% data: https:././en.wikipedia.org./wiki./Simpson%27s_rule#Composite_Simpson's_rule_for_irregularly_spaced_data
%
% © Rushikesh Kamalapurkar
%
function weights = genSimpsonsRuleWeights(t,dim)
    if nargin==1
        if min(size(t))>1
            error('Dimension input required if t is a matrix.')
        elseif size(t,1)~=1
            dim=1;
        else
            dim=2;
        end
    end
    Lengths = size(t,dim) - sum(isnan(t),dim);
    if min(mod(Lengths,2)) == 0
        error('All trajectory lengths need to be odd')
    end
        
    h=diff(t,1,dim);
    N=size(t,dim) - 1;
    weights = zeros(size(t));
    for i = 1:floor(N./2)
        j = 2.*i-1;
        if dim==1
            alpha_i = (2.*h(j+1,:).^3-h(j,:).^3+3.*h(j,:).*h(j+1,:).^2)./...
                (6.*h(j+1,:).*(h(j+1,:)+h(j,:)));
            alpha_i(isnan(h(j+1,:)))=0;
            beta_i = (h(j+1,:).^3+h(j,:).^3+3.*h(j+1,:).*h(j,:).*(h(j+1,:)+h(j,:)))./...
                (6.*h(j+1,:).*h(j,:));
            beta_i(isnan(h(j+1,:)))=0;
            eta_i = (2.*h(j,:).^3-h(j+1,:).^3+3.*h(j+1,:).*h(j,:).^2) ./...
                (6.*h(j,:).*(h(j+1,:)+h(j,:)));
            eta_i(isnan(h(j+1,:)))=0;
            weights(j,:) = weights(j,:)+eta_i;
            weights(j+1,:) = weights(j+1,:)+beta_i;
            weights(j+2,:) = weights(j+2,:)+alpha_i;
        else
            alpha_i = (2.*h(:,j+1).^3-h(:,j).^3+3.*h(:,j).*h(:,j+1).^2)./...
                (6.*h(:,j+1).*(h(:,j+1)+h(:,j)));
            alpha_i(isnan(h(j+1,:)))=0;
            beta_i = (h(:,j+1).^3+h(:,j).^3+3.*h(:,j+1).*h(:,j).*(h(:,j+1)+h(:,j)))./...
                (6.*h(:,j+1).*h(:,j));
            beta_i(isnan(h(j+1,:)))=0;
            eta_i = (2.*h(:,j).^3-h(:,j+1).^3+3.*h(:,j+1).*h(:,j).^2) ./...
                (6.*h(:,j).*(h(:,j+1)+h(:,j)));
            eta_i(isnan(h(j+1,:)))=0;
            weights(:,j) = weights(:,j)+eta_i;
            weights(:,j+1) = weights(:,j+1)+beta_i;
            weights(:,j+2) = weights(:,j+2)+alpha_i;
        end
    end
end