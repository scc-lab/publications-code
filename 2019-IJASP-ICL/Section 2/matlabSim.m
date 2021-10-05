function matlabSim
mkdir data
% monte carlo parameters
num = 200;                          % number of monte carlo sims
numTries = 3;                       % number of times to try to get a completed sim before regenerating new gains
Khigh = 15;                         % Upper bound on feedback gain K
Klow = 0.1;                         % Lower bound on feedback gain K
kCLhigh = 0.2;                      % Upper bound on concurrent learning gain kCL
kCLlow = 0.002;                     % Lower bound on concurrent learning gain kCL
GammaHigh = 3;                      % Upper bound on adaptation gain Gamma
GammaLow = 0.03;                    % Lower bound on adaptation gain Gamma
delThigh = 1;                       % Upper bound on integration/smoothing window delT
delTlow = 0.01;                     % Lower bound on integration/smoothing window delT

% Constant Parameters
N = 20;                             % history stack size
stepSize = 4e-4;                    % Simulink step size
noiseSTDev = 0.3;                   % standard deviation of AWG noise
tspan = [0:stepSize:100]';          % Simulation time span

% initial conditions
x0 = zeros(2,1);
thetaHat0 = zeros(4,1);
scriptY0 = zeros(2,4);
scriptU0 = zeros(2,1);
states0 = [x0; thetaHat0; scriptY0(:); scriptU0];

%% Sim
parfor ii = 1:num
    t1 = tic;
    doneInt = false;
    doneDeriv = false;
    while ~(doneInt && doneDeriv)
        % Gains
        K = (rand(1)*(Khigh-Klow)+Klow)*eye(2);                         % feedback gain
        kCL = (rand(1)*(kCLhigh-kCLlow)+kCLlow);                        % concurrent learning gain
        Gamma = (rand(1)*(GammaHigh-GammaLow)+GammaLow)*eye(4);         % Adaptation gain

        delT = (rand(1)*(delThigh-delTlow)+delTlow);                    % integration/smoothing window size
        buffSize = round(delT/stepSize);                                % buffer size for keeping data in memory (needed for delay/smoothing)

        % Smoothing parameters
        smoothingSpan = min(buffSize,1250);
        smoothingMat = gen_smoothing_mat(smoothingSpan,buffSize);

        % initialize buffers
        xStack = zeros(2,N);
        xdeltStack = zeros(2,N);
        scriptUstack = zeros(2,N);
        scriptYstack = zeros(2,4,N);
        xDelayBuffer = zeros(2,buffSize);
        uDelayBuffer = zeros(2,buffSize);
        YdelayBuffer = zeros(2,4,buffSize);
        persistStruct0 = struct('xStack',xStack,'xdeltStack',xdeltStack,'scriptUstack',scriptUstack,'scriptYstack',scriptYstack,'xDelayBuffer',xDelayBuffer,'uDelayBuffer',uDelayBuffer,'YdelayBuffer',YdelayBuffer);

        % simulate
        for jj = 1:numTries
            intCL = true;                       % Perform integral CL
            [statesInt,extraSignalsInt,statusInt] = customODEsolver(@dynamics,tspan,states0,persistStruct0,K,kCL,Gamma,delT,noiseSTDev,N,smoothingMat,stepSize,buffSize,intCL);
            if statusInt
                doneInt = true;
                break
            else
                doneInt = false;
            end
        end
        
        for jj = 1:numTries
            intCL = false;                      % Perform state derivative
            [statesDeriv,extraSignalsDeriv,statusDeriv] = customODEsolver(@dynamics,tspan,states0,persistStruct0,K,kCL,Gamma,delT,noiseSTDev,N,smoothingMat,stepSize,buffSize,intCL);
            if statusDeriv
                doneDeriv = true;
                break
            else
                doneDeriv = false;
            end
        end
        
        if (doneInt && doneDeriv)
            % Unpack states
            xInt = statesInt(:,1:2);
            thetaHatInt = statesInt(:,3:6);
            currSVDInt = cell2mat(extraSignalsInt(2:end));
            xDeriv = statesDeriv(:,1:2);
            thetaHatDeriv = statesDeriv(:,3:6);
            currSVDDeriv = cell2mat(extraSignalsDeriv(2:end));

            % Save data
            params = struct('K',K,'kCL',kCL,'Gamma',Gamma,'delT',delT,'N',N,'noiseSTDev',noiseSTDev);
            saveData = struct('params',params,'tspan',tspan,'xInt',xInt,'thetaHatInt',thetaHatInt,'currSVDInt',currSVDInt,'xDeriv',xDeriv,'thetaHatDeriv',thetaHatDeriv,'currSVDDeriv',currSVDDeriv);
            parforSave(ii,saveData);
        end

        % plot
        % plotData(tspan,x,thetaHat,currSVD)
    end
    toc(t1);
end

%% Generate compressed data

reducedData(num) = struct('params',[],'tspan',[],'xInt',[],'thetaHatInt',[],'currSVDInt',[],'xDeriv',[],'thetaHatDeriv',[],'currSVDDeriv',[]);
skipLen = 20; % Get every Nth data point
h = waitbar(0,'Progress...');

for ii = 1:num
    waitbar(ii/(num),h)
    data = load(['data\simout',num2str(ii),'.mat']);
    saveData = data.saveData;
    reducedData(ii).params = saveData.params;
    reducedData(ii).tspan = saveData.tspan(1:skipLen:end);
    reducedData(ii).xInt = saveData.xInt(1:skipLen:end,:);
    reducedData(ii).thetaHatInt = saveData.thetaHatInt(1:skipLen:end,:);
    reducedData(ii).currSVDInt = saveData.currSVDInt(1:skipLen:end,:);
    reducedData(ii).xDeriv = saveData.xDeriv(1:skipLen:end,:);
    reducedData(ii).thetaHatDeriv = saveData.thetaHatDeriv(1:skipLen:end,:);
    reducedData(ii).currSVDDeriv = saveData.currSVDDeriv(1:skipLen:end,:);
end

close(h);
save('reducedData.mat','reducedData');

%% Analyze

ssTime = find(reducedData(1).tspan > 40,1,'first');
eRMSint = zeros(num,2);
thetaRMSint = zeros(num,4);
eRMSderiv = zeros(num,2);
thetaRMSderiv = zeros(num,4);

% nominal values
[xd,~] = desired_trajectory(reducedData(1).tspan);
[theta,~,~] = dyn(0,[0,0]',[0,0]');

h = waitbar(0,'Progress...');
for ii = 1:num
    waitbar(ii/(num),h)
    
    % integral CL
    error = reducedData(ii).xInt - xd;
    thetaError = bsxfun(@minus,theta',reducedData(ii).thetaHatInt);
    eRMSint(ii,:) = sqrt(mean(error(ssTime:end,:).^2,1));
    thetaRMSint(ii,:) = sqrt(mean(thetaError(ssTime:end,:).^2,1));
    
    % derivative CL
    error = reducedData(ii).xDeriv - xd;
    thetaError = bsxfun(@minus,theta',reducedData(ii).thetaHatDeriv);
    eRMSderiv(ii,:) = sqrt(mean(error(ssTime:end,:).^2,1));
    thetaRMSderiv(ii,:) = sqrt(mean(thetaError(ssTime:end,:).^2,1));
    
end
close(h);

mean(eRMSint,1)
mean(thetaRMSint,1)
[minErrorInt,ind] = min(sum(eRMSint,2),[],1)
[minThetaErrorInt,ind] = min(sum(thetaRMSint,2),[],1)

goodData = ~(any(isnan(eRMSderiv),2) | any(isnan(thetaRMSderiv),2));
numGoodRuns = sum(goodData)
mean(eRMSderiv(goodData,:),1)
mean(thetaRMSderiv(goodData,:),1)
[minErrorDeriv,ind] = min(sum(eRMSderiv(goodData,:),2),[],1)
[minThetaErrorDeriv,ind] = min(sum(thetaRMSderiv(goodData,:),2),[],1)

end

function [statesDot,persistStruct,extraSignals] = dynamics(t,states,persistStruct,K,kCL,Gamma,delT,noiseSTDev,N,smoothingMat,stepSize,buffSize,intCL)

% Unpack states
x = states(1:2);
thetaHat = states(3:6);
scriptY = reshape(states(7:14),2,4);
scriptU = states(15:16);
xNoisy = x + noiseSTDev*randn(2,1);

% Unpack concurrent learning history stack
xStack = persistStruct.xStack;
xdeltStack = persistStruct.xdeltStack;
scriptUstack = persistStruct.scriptUstack;
scriptYstack = persistStruct.scriptYstack;
xDelayBuffer = persistStruct.xDelayBuffer;
uDelayBuffer = persistStruct.uDelayBuffer;
YdelayBuffer = persistStruct.YdelayBuffer;

% Auxiliary signals
[xd,xdDot] = desired_trajectory(t);
error = xNoisy - xd;

% Control
[~,Y,~] = dyn(t,xNoisy,zeros(2,1));
u = xdDot - Y*thetaHat - K*error;

% Fill buffers
xDelayBuffer = [xDelayBuffer(:,2:end), xNoisy];
uDelayBuffer = [uDelayBuffer(:,2:end), u];
YdelayBuffer = cat(3,YdelayBuffer(:,:,2:end),Y);

% State derivative estimation
if (~intCL) && (kCL > 0) && (t > delT)
    index = round(buffSize/2);
    xSmooth = smoothingMat*xDelayBuffer';
    delXsmooth = xSmooth(end,:) - xSmooth(1,:);
    xdot = (delXsmooth./(2*stepSize))';
    scriptY = YdelayBuffer(:,:,index); % careful here, scriptY used for code reuse
    scriptU = uDelayBuffer(:,index);
else
    xdot = [0;0];
end

% Singular value maximization (history stack augmenting and pruning)
if t > delT && kCL > 0
    % get current minimum singular value
    temp1 = scriptYstack;
    temp1 = reshape(permute(temp1,[2,1,3]),4,[],1)';
    minsvd = min(svd(temp1'*temp1));
    svdList = zeros(N,1);
    
    % determine singular value for each replaced data point
    for ii = 1:N
        temp2 = temp1;
        temp2(2*ii-1:2*ii,:) = scriptY;
        svdList(ii) = min(svd(temp2'*temp2));
    end
    
    % Replace data if singular value is increased
    [newsvd,ind] = max(svdList);
    if newsvd > minsvd
        if intCL
            xStack(:,ind) = xNoisy;
            xdeltStack(:,ind) = xDelayBuffer(:,1);
        else
            xStack(:,ind) = xdot;
            xdeltStack(:,ind) = zeros(2,1);
        end
        scriptUstack(:,ind) = scriptU;
        scriptYstack(:,:,ind) = scriptY;
        currSvd = newsvd;
    else
        currSvd = minsvd;
    end
else
    currSvd = 0;
end

% State derivatives
[~,~,xDot] = dyn(t,x,u);
terms = zeros(4,1,N);
for ii = 1:N
    terms(:,:,ii) = scriptYstack(:,:,ii)'*(xStack(:,ii) - xdeltStack(:,ii) - scriptUstack(:,ii) - scriptYstack(:,:,ii)*thetaHat);
end
% temp1 = mtimesx(scriptYstack,thetaHat);
% temp2 = xStack - xdeltStack - scriptUstack - temp1;
% terms = mtimesx(scriptYstack,'T',temp2);
thetaHatDot = Gamma*Y'*error + kCL*Gamma*sum(terms,3);
scriptUdot = u - uDelayBuffer(:,1);
scriptYdot = Y - YdelayBuffer(:,:,1);
statesDot = [xDot; thetaHatDot; scriptYdot(:); scriptUdot];

% repack persistent structure
persistStruct.xStack = xStack;
persistStruct.xdeltStack = xdeltStack;
persistStruct.scriptUstack = scriptUstack;
persistStruct.scriptYstack = scriptYstack;
persistStruct.xDelayBuffer = xDelayBuffer;
persistStruct.uDelayBuffer = uDelayBuffer;
persistStruct.YdelayBuffer = YdelayBuffer;

% extra signals to keep
extraSignals = currSvd;

end

function [xd,xdDot] = desired_trajectory(t)

lambda = 0.1;
if length(t) > 1
    xd = 10*bsxfun(@times,(1-exp(-lambda*t)),[sin(2*t), 0.4*cos(3*t)]);
    xdDot = 10*bsxfun(@times,(1-exp(-lambda*t)),[2*cos(2*t), -1.2*sin(3*t)])+30*bsxfun(@times,(lambda*exp(-lambda*t)),[sin(2*t), 0.4*cos(3*t)]);
else
    xd = 10*(1-exp(-lambda*t))*[sin(2*t), 0.4*cos(3*t)]';
    xdDot = 10*(1-exp(-lambda*t))*[2*cos(2*t), -1.2*sin(3*t)]'+30*(lambda*exp(-lambda*t))*[sin(2*t), 0.4*cos(3*t)]';
end

end

function smoothingMat = gen_smoothing_mat(smoothingSpan,len)

smoothingMat = zeros(len,len);
sat = @(x, lb, ub) min(max(x, lb), ub);
for ii = 1:len
    leftSpaces = ii - sat(ii-floor(smoothingSpan/2),1,len);
    rightSpaces = sat(ii+floor(smoothingSpan/2),1,len) - ii;
    spaces = min(leftSpaces,rightSpaces);
    smoothingMat(ii,ii-spaces:ii+spaces) = 1/(2*spaces+1);
end

% Reduce smoothing matrix to just relevant parts
index = round(len/2);
smoothingMat = smoothingMat(index-1:index+1,:);
smoothingMat = sparse(smoothingMat);

end

function [states,extraSignals,status] = customODEsolver(func,tspan,states0,persistStruct0,varargin)

lenT = length(tspan);
states = zeros(lenT,length(states0));
states(1,:) = states0';
extraSignals = cell(lenT,1);
persistStruct = persistStruct0;

try
    for ii = 2:lenT
        [statesDot,persistStruct,extraSignals{ii}] = func(tspan(ii-1),states(ii-1,:)',persistStruct,varargin{:});
        states(ii,:) = states(ii-1,:) + statesDot'*(tspan(ii)-tspan(ii-1));
    end
    assert(~(any(any(isnan(states))) || any(any(isinf(states(ii,:))))));
    status = true;
catch
    status = false;
%     % Unpack states
%     x = states(1:ii-1,1:2);
%     thetaHat = states(1:ii-1,3:6);
%     scriptY = reshape(states(1:ii-1,7:14)',2,4,[]);
%     scriptU = states(1:ii-1,15:16);
%     currSVD = cell2mat(extraSignals(2:ii-1));
% 
%     plotData(tspan,x,thetaHat,currSVD)
end

end

function plotData(tspan,x,thetaHat,currSVD)

% Auxiliary signals
xd = desired_trajectory(tspan);
[theta,~,~] = dyn(0,zeros(2,1),zeros(2,1));
error = x - xd;
thetaError = bsxfun(@minus,theta',thetaHat);
% scriptY = reshape(states(:,7:14)',2,4,[]);
% scriptU = states(:,15:16);

% Plot results
figure(1)
plot(tspan,error)

figure(2)
plot(tspan,thetaError)

figure(3)
plot(tspan(2:end),currSVD);

end

function [theta,Y,xDot] = dyn(t,x,u)

theta = 5*[1,2,3,4]';
Y = [x(1)^2, sin(x(2)), 0, 0;
    0, x(2)*sin(t), x(1), x(2)*x(1)];
xDot = Y*theta + u;

end