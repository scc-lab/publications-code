%% Parameters

ScalingFactor = 1;
PlotModes = 0; % Show Liouville modes
PrintModes = 0;% Save Liouville modes
TotalModesToDisplay = 5; % Number of modes to display
sortmodes='custom'; % coefficients, eigenvalues, norm, or normcoeff
ShowVideo = 1; % Show reconstruction animation
SaveVideo = 0; % Save reconstruction animation
PrintReconstruction = 0; % Save reconstruction plots

%% Explicit kernel functions that can take 3D arrays as inputs (R2020b or
% newer version of MATLAB needed)
WhichKernel='Exponential';
if isequal(WhichKernel,'Gaussian')
    mu = 500;
    K = @(X,Y) exp(-1/mu*(pagetranspose(sum(X.^2,1)) + sum(Y.^2,1) - ...
        2*pagemtimes(X,'transpose',Y,'none')));
elseif isequal(WhichKernel,'Exponential')
    mu = 1000;
    K = @(X,Y) exp(1/mu*pagemtimes(X,'transpose',Y,'none'));
elseif isequal(WhichKernel,'Linear')
    mu = 1; % Needs tuning, does not work
    K = @(X,Y) 1/mu*pagemtimes(X,'transpose',Y,'none');
end

%% Load and format data for DMD
% Trajectories can be of different lengths and irregularly sampled, but the
% lengths have to be odd.
% Inputs:
%          1) Trajectories dataset
%          2) Sample time matrix (sample times in columns)
% Trajectories Dataset Format:
%          First dimension: State 
%          Second dimension: Time (size = length of longest trajectory)
%          Third dimension: Trajectory number
%%%% Shorter trajectories need to be padded with zeros. %%%%%
% Sample Time Matrix Format:
%          First dimension: Time (size = length of longest trajectory)
%          Second dimension: Trajectory number
%%%% Sample times for shorter trajectories need to be padded with NaNs %%%%

% Load data
% Download 
% http://dmdbook.com/CODE.zip and http://dmdbook.com/DATA.zip
% and extract their contents to the same directory as this script and the
% function LiouvilleDMD.m

load('DATA/FLUIDS/CYLINDER_ALL.mat');
h = 0.02; % Sample time from DMD book
Width = 449; % Width from DMD book
Snapshots = zeros(1,89351,151);
Snapshots(1,:,:) = VORTALL;
Dimension = length(Snapshots(1,:,1));
TotalLength = length(Snapshots(1,1,:));

% Generate trajectories dataset and sample time matrix
DataSet = 5; % either 'mix' or equal to desired trajectory length (odd)
if isequal(DataSet,'mix')
    TrajectoryLength1 = 5;
    TrajectoryLength2 = 3;
    TotalTrajectories1 = 50;
    TotalTrajectories2 = 50;
    MaxTrajectoryLength = max(TrajectoryLength1,TrajectoryLength2);

    Trajectories1 = zeros(Dimension,MaxTrajectoryLength,TotalTrajectories1);
    SampleTime1 = NaN*ones(MaxTrajectoryLength,TotalTrajectories1);
    for j = 1:TotalTrajectories1
        Trajectories1(:,1:TrajectoryLength1,j) = Snapshots(1,:,j:j+(TrajectoryLength1-1));
        SampleTime1(1:TrajectoryLength1,j) = h*(0:1:TrajectoryLength1-1);
    end
    Trajectories2 = zeros(Dimension,MaxTrajectoryLength,TotalTrajectories2);
    SampleTime2 = NaN*ones(MaxTrajectoryLength,TotalTrajectories2);
    StartSecondSet = 60;
    for j = 1:TotalTrajectories2
        Trajectories2(:,1:TrajectoryLength2,j) =...
            Snapshots(1,:,j+StartSecondSet:j+StartSecondSet+(TrajectoryLength2-1));
        SampleTime2(1:TrajectoryLength2,j) =...
            h*(0:1:TrajectoryLength2-1);
    end
    Trajectories = cat(3,Trajectories1,Trajectories2);
    SampleTime = cat(2,SampleTime1,SampleTime2);
else
    TrajectoryLength = DataSet;
    TotalTrajectories = TotalLength - TrajectoryLength + 1;
    Trajectories = zeros(Dimension,TrajectoryLength,TotalTrajectories);
    SampleTime = NaN*ones(TrajectoryLength,TotalTrajectories);
    for j = 1:TotalTrajectories
        Trajectories(:,1:TrajectoryLength,j) = Snapshots(1,:,j:j+(TrajectoryLength-1));
        SampleTime(1:TrajectoryLength,j) = h*(0:1:TrajectoryLength-1);
    end
end
%% Liouville DMD
[Modes,Eigenvalues,EigenfunctionCoeffs,ReconstructionFunction] = LiouvilleDMD(K,Trajectories,SampleTime,ScalingFactor);

%% Video
if ShowVideo || SaveVideo
    if SaveVideo
        filename = ['liouville-scale-' num2str(ScalingFactor) '-reconstruction.mp4'];
        v = VideoWriter(filename, 'MPEG-4');
        v.FrameRate = 10;
        open(v);
    end
    x = 1:1:Width;
    y = 1:1:Dimension/Width;
    [XX,YY] = meshgrid(x,y);
    fig1 = figure('units','pixels','position',[100 100 1300 1000]);
    load('CH02_FLUIDS/CCcool.mat')
    colormap(CC);
    InitialState = Trajectories(1:end,1,1);
    %sgtitle('Liouville DMD reconstruction','FontSize',20);
    for i = 1:TotalLength
        t = (i-1)*h;
        ZZ = reshape(real(ReconstructionFunction(t,InitialState)),[],Width);
        subplot(2,1,1)
        fig1 = surf(XX,YY,real(ZZ));
        set(fig1,'edgecolor','none');
        set(gca,'position',[0 0.5 1 0.5]);
        view(0,90);
        colorbar
        title(['Reconstructed  $x(' num2str(t) ')$'],'FontSize',20,'Interpreter','latex','Units', 'normalized', 'Position', [0.5, 0.9, 0]);
        set(gca,'YTickLabel',[],'XTickLabel',[],'XColor', 'none','YColor','none');
        subplot(2,1,2)
        ZZ = reshape(Snapshots(1,:,i),[],Width);
        fig1 = surf(XX,YY,real(ZZ));
        set(fig1,'edgecolor','none');
        set(gca,'position',[0 0 1 0.5]);
        view(0,90);
        colorbar
        title(['True  $x(' num2str(t) ')$'],'FontSize',20,'Interpreter','latex','Units', 'normalized', 'Position', [0.5, 0.9, 0]);
        set(gca,'YTickLabel',[],'XTickLabel',[],'XColor', 'none','YColor','none');
        drawnow;
        if SaveVideo
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
    end
    if SaveVideo
        close(v);
    end
end
%% Print reconstruction
if PrintReconstruction
    x = 1:1:Width;
    y = 1:1:Dimension/Width;
    [XX,YY] = meshgrid(x,y);
    ReconstructionIndices = [10 55 80 113 150];
    load('CH02_FLUIDS/CCcool.mat') 
    InitialState = Trajectories(1:end,1,1);
    for ii=1:numel(ReconstructionIndices)
        i = ReconstructionIndices(ii);
        t = (i-1)*h;
        ZZ = reshape(real(ReconstructionFunction(t,InitialState)),[],Width);
        fig1 = surf(XX,YY,real(ZZ));
        set(fig1,'edgecolor','none');
        view(0,90);
        colormap(CC); 
        set(gca,'position',[0 0 1 1],'YTickLabel',[],'XTickLabel',[],'XColor', 'none','YColor','none');
        if(PrintReconstruction == 1)
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperSize', [6 5]);
            set(gcf, 'PaperPosition', [0 0 6 5]);
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
            filename = ['liouville-scale-' num2str(ScalingFactor) '-reconstruction-at-' num2str(t) '.pdf'];
            saveas(gcf,filename);
        end
        ZZ = reshape(Snapshots(1,:,i),[],Width);
        fig1 = surf(XX,YY,real(ZZ));
        set(fig1,'edgecolor','none');
        view(0,90);
        set(gca,'position',[0 0 1 1],'YTickLabel',[],'XTickLabel',[],'XColor', 'none','YColor','none');
        axis tight
        if(PrintReconstruction == 1)
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperSize', [6 5]);
            set(gcf, 'PaperPosition', [0 0 6 5]);
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
            filename = ['liouville-scale-' num2str(ScalingFactor) '-true-at-' num2str(t) '.pdf'];
            saveas(gcf,filename);
        end
    end
end

%% Plot modes
if PlotModes
    TotalTrajectories = size(Trajectories,3);
    load('CH02_FLUIDS/CCcool.mat') 
    TotalModesToDisplay = min(TotalModesToDisplay,TotalTrajectories);
    if isequal(sortmodes,'coefficients')
        [~,sortindex] = sort(abs(EigenfunctionCoeffs),'descend');
    elseif isequal(sortmodes,'eigenvalues')
        [~,sortindex] = sort(real(Eigenvalues),'descend');
    elseif isequal(sortmodes,'norm')
        [~,sortindex] = sort(vecnorm(Modes),'descend');
    elseif isequal(sortmodes,'normcoeff')
        [~,sortindex] = sort(abs(EigenfunctionCoeffs).*vecnorm(Modes),'descend');
    elseif isequal(sortmodes,'custom')
        modesofinterest=[142 120 103 94 80 65 51 36 20 41 2 18];
        sortindex=[modesofinterest setdiff(1:TotalTrajectories,modesofinterest)];
    else
        sortindex = 1:TotalTrajectories;
    end
    for i = 1:TotalModesToDisplay
        s = sortindex(i);
        Z = reshape(Modes(:,s),[],Width);
        x = 1:1:Width;
        y = 1:1:length(Z(:,1));
        [X,Y] = meshgrid(x,y);
        figure
        fig1 = surf(X,Y,real(Z));
        set(fig1,'edgecolor','none');
        view(0,90)
        colormap(CC); 
        %title("Real Part of Liouville Mode " + s,'fontsize',8);
        if(PrintModes == 1)
            set(gca,'position',[0 0 1 1],'YTickLabel',[],'XTickLabel',[],'XColor', 'none','YColor','none');
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperSize', [6 5]);
            set(gcf, 'PaperPosition', [0 0 6 5]);
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
            filename = ['rbf-liouville-scale-' num2str(ScalingFactor) '-real-' num2str(i) '-mode-' num2str(s) '.pdf'];
            saveas(gcf,filename);
        end
        
        figure
        fig2 = surf(X,Y,imag(Z));
        set(fig2,'edgecolor','none');
        view(0,90)
        colormap(CC); 
        %title("Imaginary Part of Liouville Mode " + s,'fontsize',8);
        if(PrintModes == 1)
            set(gca,'position',[0 0 1 1],'YTickLabel',[],'XTickLabel',[],'XColor', 'none','YColor','none');
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperSize', [6 5]);
            set(gcf, 'PaperPosition', [0 0 6 5]);
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
            filename = ['rbf-liouville-scale-' num2str(ScalingFactor) '-imag-' num2str(i) '-mode-' num2str(s) '.pdf'];
            saveas(gcf,filename);
        end
    end   
figHandles = findobj('Type', 'figure');
f_figureplace(figHandles,5,6);
end