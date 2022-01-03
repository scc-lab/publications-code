function [OutputVector] = evaluate_occupation_kernels(x,approximate_paths,mu,h,SimpsonsRuleVector)
% This gives a column vector of the occupation kernels corresponding to a
% collection of paths in approximat_paths evaluated at x.

TotalLength = length(approximate_paths(1,:,1));

OutputVector = zeros(length(approximate_paths(1,1,:)),1);

for ii = 1:length(approximate_paths(1,1,:))
    XY = 2/mu.*squeeze(approximate_paths(:,:,ii))'*squeeze(x);
%     XX = -1/mu.*diag(approximate_paths(:,:,ii)'*approximate_paths(:,:,ii));
%     YY = -1/mu.*(diag(diag(squeeze(x)'*squeeze(x)))*ones(TotalLength,1)); % This ones matrix might be an issue.

    OutputVector(ii) = h/3*(SimpsonsRuleVector*exp(XY));
end

end

