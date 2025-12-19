% This script generates snapshots of the position of mesh points in an
% undamped cantilever beam and uses second order LiouvilleDMD to generate a 
% predictive model from the data.
%
% © Rushikesh Kamalapurkar, Ben Russo, and Joel Rosenfeld
%
clear all
close all
addpath('../../lib');

%% Parameters
directReconstruction = 0; % Do direct reconstruction (indirect by default)
showVideo = 0; % Show reconstruction animation
saveVideo = 0; % Save animation as mp4 video
printReconstruction = 0; % Save reconstruction snapshots
printGeneralization = 0; % Save generalization and comparison plots
printError = 0; % Save RMS error plots

%% Kernel 
% (R2020b or newer version of MATLAB needed)
% Second order
mu2 = 0.0005; % 0.08 for exponential, 0.0005 for linear
K2 = KernelRKHS('Linear',mu2);
l2 = 1e-14; % Regularization parameter for rank deficient Gram matrices

% First order
mu1 = 500; % 1000 for exponential, 0.0005 for linear
K1 = KernelRKHS('Linear',mu1);
l1 = 1e-10; % Regularization parameter for rank deficient Gram matrices
% 3e-12 for exponential

%% Generate and format data for DMD
% Generate data
% See https://www.mathworks.com/help/pde/ug/dynamics-of-a-damped-cantilever-beam.html
% MATLAB Partial Differential Equations Toolbox required
[V,msh,longestPeriod] = cantileverBeamUndampedTransient(10,0);
h = longestPeriod/100;

% Segment single trajectory into multiple trajectories
OutputSnapshots = zeros(1,size(msh.Nodes,1)*size(V.Displacement.x,1),size(V.Displacement.x,2));
OutputSnapshots(1,:,:) = [V.Displacement.x;V.Displacement.y];
StateSnapshots(1,:,:) = [V.Displacement.x;V.Displacement.y;V.Velocity.x;V.Velocity.y];
Velocities = zeros(1,size(msh.Nodes(1,:),1)*size(V.Displacement.x,1),size(V.Displacement.x,2));
outputDimension = length(OutputSnapshots(1,:,1));
order = 2;
stateDimension = outputDimension * order;
TotalLength = length(OutputSnapshots(1,1,:));
NumNodes = size(msh.Nodes,2);
% Generate trajectories dataset and sample time matrix
TrajectoryLength = 31;
TotalTrajectories = TotalLength - TrajectoryLength + 1;
Output = zeros(outputDimension,TrajectoryLength,TotalTrajectories);
State = zeros(stateDimension,TrajectoryLength,TotalTrajectories);
Derivatives = zeros(outputDimension,TotalTrajectories);
SampleTime = NaN*ones(TrajectoryLength,TotalTrajectories);
for j = 1:TotalTrajectories
    Output(:,1:TrajectoryLength,j) = OutputSnapshots(1,:,j:j+(TrajectoryLength-1));
    State(:,1:TrajectoryLength,j) = StateSnapshots(1,:,j:j+(TrajectoryLength-1));
    Derivatives(:,j) = Velocities(:,j);
    SampleTime(1:TrajectoryLength,j) = h*(0:1:TrajectoryLength-1);
end

%% Liouville DMD
% Second order
[~,~,~,r2,f2] = SecondOrderLiouvilleDMD(K2,Output,SampleTime,[],l2);

% First order
[~,~,~,r1,f1] = LiouvilleDMD(K1,State,SampleTime,[],l1);

%% Parameters for plotting
ll = [0,0.006,0.011,0.019,0.0275,0.037,0.045,0.0555,0.064,0.073,0.0825,0.091,0.1];
hl = [0,0.01,0.018,0.0275,0.036,0.044,0.0551,0.063,0.0722,0.0824,0.0905,0.095,0.1];
beamWidth = numel(ll);
beamLength = 490;
[YY,II] = sort(msh.Nodes(2,:));
indexArray = zeros(1,beamLength*beamWidth);
for i=1:beamWidth
    temp = find((YY <= hl(i)) .* (YY >= ll(i)));
    indexArray((i-1)*beamLength+1:i*beamLength) = temp(1:beamLength);
end
II = II(indexArray);

%% Indirect reconstruction
f22 = @(t,x) [x(2*NumNodes+1:end); f2(x(1:2*NumNodes))];
initialState = [Output(1:end,1,1);Derivatives(:,1)];
Time = ((1:size(V.Displacement.x,2))-1)*h;
[~,zi2] = ode45(@(t,x) f22(t,x), Time, initialState);

%% Plots
if showVideo || saveVideo
    if saveVideo
        filename = ['secondOrderLiouvilleBeamIndirect.mp4'];
        v = VideoWriter(filename, 'MPEG-4');
        v.FrameRate = 10;
        open(v);
    end
    figure('units','pixels','position',[100 100 1300 1000]);
    for i = 1:150
        ey2 = zi2(i,NumNodes+1:2*NumNodes).' - V.Displacement.y(:,i);
        ex2 = zi2(i,1:NumNodes).' - V.Displacement.x(:,i);
        pointwise_error_norm = sqrt(ex2.^2 + ey2.^2);
        newNodes = zeros(size(msh.Nodes));
        newNodes(1,:) = msh.Nodes(1,:)+V.Displacement.x(:,i).';
        newNodes(2,:) = msh.Nodes(2,:)+V.Displacement.y(:,i).';
        pdeplot(newNodes,msh.Elements,"XYData",pointwise_error_norm,"ColorMap","parula")
        colorbar;
        clim([0 2.5e-3]);
        ylim([-0.02,0.12]);
        xlim([1,5]);
        drawnow;
    end
end

if printReconstruction
        ReconstructionIndices = [1 30 70 100 130];
        fig1 = figure();
        for ii=1:numel(ReconstructionIndices)
            i = ReconstructionIndices(ii);
            ey2 = zi2(i,NumNodes+1:2*NumNodes).' - V.Displacement.y(:,i);
            ex2 = zi2(i,1:NumNodes).' - V.Displacement.x(:,i);
            error_norm_pointwise = sqrt(ex2.^2 + ey2.^2);
            subplot('position',[0 0.95*(ii-1)/numel(ReconstructionIndices) 0.8 0.95/numel(ReconstructionIndices)],'parent',fig1);
            newNodes = zeros(size(msh.Nodes));
            newNodes(1,:) = msh.Nodes(1,:)+V.Displacement.x(:,i).';
            newNodes(2,:) = msh.Nodes(2,:)+V.Displacement.y(:,i).';
            pdeplot(newNodes,msh.Elements,"XYData",error_norm_pointwise,"ColorMap","parula","ColorBar","off")
            clim([0 2.5e-3]);
            ylim([-0.02,0.12]);
            xlim([1,5]);
            set(gca,'YTickLabel',[],'XTickLabel',[]);
            title({'t =', [num2str(t*1000,2) ' ms']},'FontSize',11,'fontweight','normal','Units','normalized','Position',[1.06, 0.3, 0],'Interpreter','LaTeX');
        end
        ax1 = axes(fig1,'visible','off');
        c = colorbar(ax1,'Position',[0.9 0.1 0.022 0.8],'TickLabelInterpreter','latex');
        clim([0 2.5e-3]);
        set(ax1,'fontsize',12,'TickLabelInterpreter','latex');
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [7 5]);
        set(gcf, 'PaperPosition', [0 0 7 5]);
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        fig1.Renderer='painters';
        filename = 'secondOrderLiouvilleBeamIndirectReconstruction';
        saveas(fig1,filename,'pdf');
end

if printError
    supNorm = max(norm([V.Displacement.x(:); V.Displacement.y(:)]));
    fig1 = figure();
    relativeError = zeros(size(Time));
    for i = 1:size(V.Displacement.x,2)
        ey2 = zi2(i,NumNodes+1:2*NumNodes).' - V.Displacement.y(:,i);
        ex2 = zi2(i,1:NumNodes).' - V.Displacement.x(:,i);
        error_norm = norm([ex2 ey2]);
        relativeError(i) = error_norm/supNorm;
    end
    plot(Time,relativeError,'LineWidth',2)
    ylabel('RMS Error','Interpreter','LaTeX');
    xlabel('Time [s]','Interpreter','LaTeX');
    set(gca,'fontsize',12,'TickLabelInterpreter','latex');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 5]);
    set(gcf, 'PaperPosition', [0 0 6 5]);
    filename = 'secondOrderLiouvilleBeamIndirectReconstructionError';
    saveas(fig1,filename,'pdf');
end

% Generalization
[VG,~,~] = cantileverBeamUndampedTransient(15,0);
initialState = [VG.Displacement.x(:,1);
                VG.Displacement.y(:,1);
                zeros(2*NumNodes,1)];
Time = ((1:size(VG.Displacement.x,2))-1)*h;
[~,zi1g] = ode45(@(t,x) f1(x), Time, initialState);
[~,zi2g] = ode45(@(t,x) f22(t,x), Time, initialState);
supNorm = max(norm([VG.Displacement.x(:); VG.Displacement.y(:)]));
fig1 = figure();
ei1g = zeros(size(Time));
ei2g = zeros(size(Time));
for i = 1:size(VG.Displacement.x,2)
    error_norm2 = norm([zi2g(i,1:NumNodes).' - VG.Displacement.x(:,i),...
        zi2g(i,NumNodes+1:2*NumNodes).' - VG.Displacement.y(:,i)]);
    ei2g(i) = error_norm2/supNorm;
    error_norm1 = norm([zi1g(1,1:NumNodes).' - VG.Displacement.x(:,i),...
        zi1g(1,NumNodes+1:2*NumNodes).' - VG.Displacement.y(:,i)]);
    ei1g(i) = error_norm1/supNorm;
end
temp = [Time.', ei1g.', ei2g.'];
save('LiouvilleBeamIndirectGeneralizationError.dat','temp','-ascii');
plot(Time,ei1g,Time,ei2g,'LineWidth',1.5)
ylabel('Relative RMS Error','Interpreter','LaTeX');
xlabel('Time [s]','Interpreter','LaTeX');
l=legend('First order DMD','Second order DMD');
set(l,'Interpreter','latex');
set(gca,'fontsize',12,'TickLabelInterpreter','latex');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
if printGeneralization
    filename = 'LiouvilleBeamIndirectGeneralizationError';
    saveas(fig1,filename,'pdf');
end

%% Plots using direct reconstruction
if directReconstruction
    if showVideo || saveVideo
        if saveVideo
            filename = 'secondOrderLiouvilleBeamDirect.mp4';
            v = VideoWriter(filename, 'MPEG-4');
            v.FrameRate = 10;
            open(v);
        end
        figure('units','pixels','position',[100 100 1300 1000]);
        for i = 1:150
            t = (i-1)*h;
            zd2 = r2(t,Output(1:end,1,1),Derivatives(:,1));
            ey2 = zd2(NumNodes+1:end) - V.Displacement.y(:,i);
            ex2 = zd2(1:NumNodes) - V.Displacement.x(:,i);
            error_norm_pointwise = sqrt(ex2.^2 + ey2.^2);
            newNodes = zeros(size(msh.Nodes));
            newNodes(1,:) = msh.Nodes(1,:)+V.Displacement.x(:,i).';
            newNodes(2,:) = msh.Nodes(2,:)+V.Displacement.y(:,i).';
            pdeplot(newNodes,msh.Elements,"XYData",error_norm_pointwise,"ColorMap","parula")
            colorbar;
            clim([0 2.5e-3]);
            ylim([-0.02,0.12]);
            xlim([1,5]);
            drawnow;
            if saveVideo
                frame = getframe(gcf);
                writeVideo(v,frame);
            end
        end
        if saveVideo
            close(v);
        end
    end

    if printReconstruction
        ReconstructionIndices = [1 30 70 100 130];
        fig1 = figure();
        for ii=1:numel(ReconstructionIndices)
            i = ReconstructionIndices(ii);
            t = (i-1)*h;
            zd2 = r2(t,Output(1:end,1,1),Derivatives(:,1));
            ey2 = zd2(NumNodes+1:end) - V.Displacement.y(:,i);
            ex2 = zd2(1:NumNodes) - V.Displacement.x(:,i);
            error_norm_pointwise = sqrt(ex2.^2 + ey2.^2);
            subplot('position',[0 0.95*(ii-1)/numel(ReconstructionIndices) 0.8 0.95/numel(ReconstructionIndices)],'parent',fig1);
            newNodes = zeros(size(msh.Nodes));
            newNodes(1,:) = msh.Nodes(1,:)+V.Displacement.x(:,i).';
            newNodes(2,:) = msh.Nodes(2,:)+V.Displacement.y(:,i).';
            pdeplot(newNodes,msh.Elements,"XYData",error_norm_pointwise,"ColorMap","parula","ColorBar","off")
            clim([0 2.5e-3]);
            ylim([-0.02,0.12]);
            xlim([1,5]);
            set(gca,'YTickLabel',[],'XTickLabel',[]);
            title({'t =', [num2str(t*1000,2) ' ms']},'FontSize',11,'fontweight','normal','Units','normalized','Position',[1.06, 0.3, 0],'Interpreter','LaTeX');
        end
        ax1 = axes(fig1,'visible','off');
        c = colorbar(ax1,'Position',[0.9 0.1 0.022 0.8],'TickLabelInterpreter','latex');
        clim([0 2.5e-3]);
        set(ax1,'fontsize',12,'TickLabelInterpreter','latex');
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [7 5]);
        set(gcf, 'PaperPosition', [0 0 7 5]);
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        fig1.Renderer='painters';
        filename = 'secondOrderLiouvilleBeamDirectReconstruction';
        saveas(fig1,filename,'pdf');
    end

    if printError
        supNorm = max(norm([V.Displacement.x(:); V.Displacement.y(:)]));
        fig1 = figure();
        Time = ((1:size(V.Displacement.x,2))-1)*h;
        relativeError = zeros(size(Time));
        for i = 1:size(V.Displacement.x,2)
            t = (i-1)*h;
            zd2 = r2(t,Output(1:end,1,1),Derivatives(:,1));
            ey2 = zd2(NumNodes+1:end) - V.Displacement.y(:,i);
            ex2 = zd2(1:NumNodes) - V.Displacement.x(:,i);
            error_norm = norm([ex2 ey2]);
            relativeError(i) = error_norm/supNorm;
        end
        plot(Time,relativeError,'LineWidth',2)
        ylabel('RMS Error','Interpreter','LaTeX');
        xlabel('Time [s]','Interpreter','LaTeX');
        set(gca,'fontsize',12,'TickLabelInterpreter','latex');
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [6 5]);
        set(gcf, 'PaperPosition', [0 0 6 5]);
        filename = 'secondOrderLiouvilleBeamDirectReconstructionError';
        saveas(fig1,filename,'pdf');
    end
end
