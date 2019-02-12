function [xH_iDot,thetaH_iDot,SIDConLear_i,SIDRank_i] = ...
    ADPCLNNTSSysID(i,xH_i,thetaH_i,u_i,t,xi,k_i,...
    GammaTh_i,kTheta_i,M_theta_i,fs,F,g)

% This simulink block implements concurrent-learning based system
% identification using a Savitzky–Golay smoothing filter for state
% derivative computation. 
% fs: Sampling frequency
% M: number of stored data points
% F: Savitzky–Golay filter window length
% g: Savitzky–Golay filter differentiation matrix

global XX_iC; % Data store memory for storing the values of x
global U_iC; % Data store memory for storing the values of u
global SIGMATHETAS_iC; % Data store for storing the values of sigth,
           % made of sigth matrices stacked in a row
global xDmGUS_iC; % Data store memory for storing the values of xDot-G*u

XX_i=XX_iC{i};
U_i=U_iC{i};
SIGMATHETAS_i=SIGMATHETAS_iC{i};
xDmGUS_i=xDmGUS_iC{i};
%Dynamics
G=cos(2*xi)+2;

% if model~=1
ite=round(t*fs)+1; % Current iteration number
xT = xi-xH_i;

% Basis functions
[sigth_i,~]=ADPCLNNTSsigmatheta(i,xi,t);

% At every time ite store the past F values of the state. Note tH 
% Q*g(:,2)/(1/fs) gives the state derivative at ite-(F-1)/2
if ite<=F
    XX_i(:,ite)=xi;
else
    XX_i = [XX_i(:,2:F) xi];
end

% At every time ite store the past (F+1)/2 values of the control. Actually, 
% we only need the value of u at ite-(F-1)/2. However, I could not find a 
% way to store just tH single value. Also, even if U is initialized as a 
% row vector, MATLAB seems to think its a column.
if ite<=(F+1)/2
    U_i(ite)=u_i;
else
    U_i = [U_i(2:(F+1)/2) u_i];
end

% The state derivative information is only available if ite>=F use
% singular value maximization algorithm to replace elements of SIGMATHETAS
% and XDmGUS with new data. At each time ite, the value of the 
% state at ite-(F-1)/2 can be found in XX(:,(F+1)/2)
index=0;
if ite>=F
    xj=XX_i(:,(F+1)/2);
    [sigth_ij,~]=ADPCLNNTSsigmatheta(i,xj,t);
    REGste1 = SIGMATHETAS_i*SIGMATHETAS_i'; % Compute Summation
    svd1=min(svd(REGste1));
    for j=1:M_theta_i
         % Store existing data in a temporary variable
        SIGMATHETASte=SIGMATHETAS_i;
        % Replac i'th sigth matrix with sigthj
        SIGMATHETASte(:,j)=sigth_ij;
        REGste2 = SIGMATHETASte*SIGMATHETASte'; % Compute Summation
        svd2=min(svd(REGste2)); % Minimum singular value after replacement
        if svd1<svd2 % If replacement results in a higher singular value,
            index=j; % note the index of replacement and
            svd1=svd2; % keep searching for an even higher singular value.           
        end
    end
    
% At this point, index points to the matrix sigth in the stack 
% SIGMATHETAS which, when replaced with the new matrix sigthj,
% results in the highest singular value.
% If index=0, the new data point is discarded. Otherwise, the stack
% SIGMATHETAS, and correspondingly the stack xDmGUS is updated with the 
% new data point.
% Note tH since there is a delay of (F-1)/2 for derivative computation,
% the "new" data point we are talking about is the one tH was recorded 
% (F-1)/2 steps ago.
    if index ~= 0
        SIGMATHETAS_i(:,index)=sigth_ij;
        uj=U_i(1,1);
        xDj=XX_i*g(:,2)/(1/fs);
        Gj = cos(2*xj)+2; 
        xDmGUS_i(:,index)=xDj-Gj*uj;
    end
end

xH_iDot = thetaH_i'*sigth_i + G*u_i + k_i*xT;

% Final computation for the concurrent learning-based component of the
% adaptive update law. The concurrent learning component is zero for the
% first F timesteps.
REGs = SIGMATHETAS_i*SIGMATHETAS_i'*(ite>=F);
SIDRank_i = min(svd(REGs)); % For plots
REG1s = SIGMATHETAS_i*xDmGUS_i'*(ite>=F);
SIDConLear_i = (REG1s-REGs*thetaH_i)/M_theta_i;
thetaH_iDot=(GammaTh_i*sigth_i*xT'+kTheta_i*GammaTh_i*SIDConLear_i);
XX_iC{i}=XX_i;
U_iC{i}=U_i;
SIGMATHETAS_iC{i}=SIGMATHETAS_i;
xDmGUS_iC{i}=xDmGUS_i;
% else
% SIDConLear_i = zeros(size(thetaH_i));
% SIDRank_i = 0;
% xH_iDot = zeros(size(xH_i));    
% thetaH_iDot = zeros(size(thetaH_i));
% end

end

