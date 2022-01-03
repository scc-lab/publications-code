%% Motion Tomography Code
function [metrics,X,Y,VectorField,true_paths,approximate_paths,VectorFieldApprox] = OccMotionTomography(numIterations,flowField)
    % When the exponential dot product is used, numerical inf is easy to
    % kick out. So make sure this is large enough to damp down your data to
    % survive the summation process, but not too big to kill your data.
    mu = 50; 

% Initializing Variables
    % Simulation Parameters
        T = 1; % Timelimit
        h = 0.01; % h spacing
        speed = 1; % Constant speed for all 

% True Flow Field Function
    
    if(flowField == 1)
        % Linear combination of Gaussians with width 3.
        F = @(x) ...
            1/8*[5*exp(-2/1*norm(x-[1/4;1/4])^2) ...
            - 0.2*exp(-1/1*norm(x-[1/4;3/4])^2) ...
            + 2*exp(-1/1*norm(x-[3/4;3/4])^2) ...
            - 5*exp(-2/1*norm(x-[3/4;1/4])^2); ...
            3*exp(-1/1*norm(x-[1/4;1/4])^2) ...
            + exp(-1/1*norm(x-[1/4;3/4])^2) ...
            - 3*exp(-3/1*norm(x-[3/4;3/4])^2) ...
            + exp(-1/1*norm(x-[3/4;1/4])^2)];
    elseif(flowField==2)
        % Spiral flow field
        F = @(x) [x(2);
            -0.2*(x(1))];
    elseif(flowField==3)
        % Constant flow field
        F = @(x) [0.2;
            0.1];
    end
    
% Initial Trajectories

    % Vector of Initial Points and Directions
    
    load('initialpoints.mat'); % rinit = [ rinit(1,1), ...; rinit(2,1), ...];
    load('direction.mat'); % direction = [ direction(1,1), ...; direction(2,1), ...];
    
    TotalPoints = length(rinit(1,:));
    % Straight line paths
    anticipated_paths = zeros(2,floor(T/h),TotalPoints);
    
    for i = 1:TotalPoints
        % rk4 with constant dynamics speed*direction(i,:) and initial point
        % rinit(i,:)
        anticipated_paths(:,:,i) = rk4method( rinit(:,i), @(x) speed*direction(:,i), h, T); % Store Trajectory Here
    end
    
    approximate_paths = anticipated_paths;

    
% True Trajectories

    % For each initial point and direction use RK4 for time 0 to T
    true_paths = zeros(2,floor(T/h),TotalPoints);

    for i = 1:TotalPoints
        % rk4 with dynamics speed*direction(i,:)+F(x) and initial point
        % rinit(i,:)
        true_paths(:,:,i) = rk4method( rinit(:,i), @(x) speed*direction(:,i) + F(x), h, T); % Store Trajectory Here
    end
    
    % Parity Check for Simpson's Rule
    
    if(mod(length(true_paths(1,:,1)),2) == 0)
        true_paths = true_paths(:,1:end-1,:);
        anticipated_paths = anticipated_paths(:,1:end-1,:);
        approximate_paths = anticipated_paths;
        old_approximate_paths = zeros(size(approximate_paths));
    end
    
    LengthTrajectories = length(true_paths(1,:,1));
    
    % Zero Initial Approximation

    approx_F = @(x) [0;0];
            
    weights_x = zeros(TotalPoints,1);
    weights_y = zeros(TotalPoints,1);


%% Iteration Loop

for m = 1:numIterations
    
    
    % We are going to use the approximate_paths to make our occupation
    % kernels. approximate_paths is a three dimensional matrix, with the
    % first dimension being the dimension of the system (2 in this case),
    % the second is the time points for the trajectory (which is indexed by
    % the third entry).
    
    % Need a Simpson's rule vector.
    TotalLength = length(squeeze(approximate_paths(1,:,1)));
    
    SimpsonsRuleVector = 1:1:TotalLength;
        SimpsonsRuleVector = 2 + (1 + (-1).^SimpsonsRuleVector);
        SimpsonsRuleVector(1) = 1;
        SimpsonsRuleVector(end) = 1;
    
    % This will be used in the generation of the Gram matrix later.
    SimpsonsRuleMatrix = SimpsonsRuleVector'*SimpsonsRuleVector; 
    
    % Difference Between Anticipated and True Endpoints
    
    D = true_paths(:,end,:) - approximate_paths(:,end,:);
    
    D = squeeze(D)';
    
    % Approximation Compensation
    
        % Interaction Matrix (like a gram matrix)
            % Rows should come from new trajectories
            % Columns should come from old trajectories
            
            InteractionMatrix = zeros(TotalPoints); % 2 dimensional square matrix.
            
            for ii = 1:TotalPoints
                for jj = 1:TotalPoints
                    % Here we are going to do something weird. We are going
                    % to make some matrices by multiplying the approximate
                    % trajectory matrices. Then we exponentiate them. This
                    % is *way* faster than doing this entry by entry.
                    % Remember: New trajectories are rows here.
                    
                    % If something breaks... transpose this...
                    
                    XY = 2/mu.*squeeze(approximate_paths(:,:,ii))'*squeeze(old_approximate_paths(:,:,jj));
%                     XX = -1/mu.*(ones(TotalLength,TotalLength))*diag(diag(approximate_paths(:,:,ii)'*approximate_paths(:,:,ii)));
%                     YY = -1/mu.*(diag(diag(squeeze(old_approximate_paths(:,:,jj))'*squeeze(old_approximate_paths(:,:,jj))))*ones(TotalLength,TotalLength));
                    
                    InteractionMatrix(ii,jj) = h^2/9*sum(sum(SimpsonsRuleMatrix.*exp(XY)));
                end
            end
                    
                    
            
    
    FCompensate = InteractionMatrix*[weights_x,weights_y]; % Need to come back and fix this. This is the inner product of the old and new occupation kernel functions.
    
    
    D = D + FCompensate;
    
    % Approximation of Flow Field with Respect to Occupation Kernels
    GramMatrix = zeros(TotalPoints);
    
            % Same deal as with the interaction matrix
            
            for ii = 1:TotalPoints
                for jj = 1:TotalPoints
                    % Here we are going to do something weird. We are going
                    % to make some matrices by multiplying the approximate
                    % trajectory matrices. Then we exponentiate them. This
                    % is *way* faster than doing this entry by entry.
                    % Remember: New trajectories are rows here.
                    
                    XY = 2/mu.*squeeze(approximate_paths(:,:,ii))'*squeeze(approximate_paths(:,:,jj));
                    %XX = -1/mu.*(ones(TotalLength,TotalLength))*diag(diag(squeeze(approximate_paths(:,:,ii))'*squeeze(approximate_paths(:,:,ii))));
                    %YY = -1/mu.*(diag(diag(squeeze(approximate_paths(:,:,jj))'*squeeze(approximate_paths(:,:,jj))))*ones(TotalLength,TotalLength));
                    
                    GramMatrix(ii,jj) = h^2/9*sum(sum(SimpsonsRuleMatrix.*exp(XY)));
                end
            end
    
    weights_x = pinv(GramMatrix)*D(:,1);
    weights_y = pinv(GramMatrix)*D(:,2);
    
    % This is the evaluation of our generated approximation. The weights
    % for x and y are multiplied by a column vector of our occupation
    % kernels evaluated at x.
    approx_F = @(x) [weights_x';weights_y']*evaluate_occupation_kernels(x,approximate_paths,mu,h,SimpsonsRuleVector);
        
    % New Anticipated Trajectories
    
    old_approximate_paths = approximate_paths;
    
    for i = 1:TotalPoints
        % rk4 with dynamics speed*direction(i,:)+approx_F(x) and initial point
        % rinit(i,:)
        holdontothis =  rk4method( rinit(:,i), @(x) speed*direction(:,i)+approx_F(x), h, T); % Store Trajectory Here
        approximate_paths(:,:,i) = holdontothis(:,1:end-1);
    end
    

% End Iteration Loop
end

%% Compare True Flow Field and Approximate Flow Field

x = 0:0.05:1;
y = 0:0.05:1;

[X,Y] = meshgrid(x,y);

X = reshape(X,[],1);
Y = reshape(Y,[],1);

% True Flow Field
VectorField = zeros(2,length(X));
for i=1:length(X)
    VectorField(:,i) = F([X(i);Y(i)]);
end

% Approximate Flow Field
VectorFieldApprox = zeros(2,length(X));
for i=1:length(X)
    VectorFieldApprox(:,i) = approx_F([X(i);Y(i)]);
end

L=length(VectorField(1,:));

%lets compute the max norm difference.
Norm_diffs=[];
for i=1:L
    Norm_diffs = [Norm_diffs, norm(VectorField(:,i)-VectorFieldApprox(:,i))];
end

percenterror=[];
for i=1:L
    percenterror = [percenterror, (norm(VectorField(:,i)-VectorFieldApprox(:,i)))/norm(VectorField(:,i))];
end
percenterror = percenterror(isnan(percenterror));
percenterror = percenterror(percenterror~=Inf);
metrics = [max(percenterror), mean(percenterror), max(Norm_diffs), mean(Norm_diffs)];
end
