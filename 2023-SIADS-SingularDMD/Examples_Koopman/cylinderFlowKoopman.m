% This code was created to accompany the manuscript "Singular Dynamic Mode
% Decompositions" by Joel A. Rosenfeld and Rushikesh Kamalapurkar. This is
% an explicit implementation of the kernel DMD method in Williams et al.,
% 2015, see for details
%
% https://doi.org/10.3934/jcd.2015005

clear all
DATAPATH = '../DATA';
addpath('../Functions')

%% Parameters
% modes
TotalModesToDisplay = 10;
sortmodes='custom'; % coefficients, eigenvalues, norm, or normcoeff
PlotModes = 0;
PrintModes = 0;
% video
ShowReconVideo = 1;
SaveReconVideo = 0;
% plots
PrintRecon = 0;
%% Kernel selection
mu = 500; 
K = Kernel('Exponential',mu); 

%% Load and format data for DMD
load([DATAPATH '/FLUIDS/CYLINDER_ALL.mat']);
DATA = VORTALL; % Using vorticity data to compare with the DMD book
Length = 150; % Number of snapshots used for DMD. Should be between 1 and
             % 150, cannot be 150 if PlotPrediction is ON.
h = 0.02; % Time step
Width = 449; % Width needed for plotting only
l = 1e-10; % Gram matrix regularization parameter

Dimension = size(DATA,1);
% Input values X
X = DATA(:,1:Length);
% Output values Y = F(X)
Y = DATA(:,2:Length+1);

%% Kernel DMD
tic
[~,~,~,ContinuousReconFun,DiscreteReconFun] = WilliamsKDMD(X,Y,K,h,l);
toc
%% Video
if ShowReconVideo || SaveReconVideo
    if SaveReconVideo
        v = VideoWriter('KDMDRecon.mp4', 'MPEG-4');
        v.FrameRate = 10;
        open(v);
    end
    x0 = X(:,1);
    x = 1:1:Width;
    y = 1:1:Dimension/Width;
    [XX,YY] = meshgrid(x,y);
    fig1 = figure('units','pixels','position',[100 50 1300 450]);
    load([DATAPATH '/CCcool.mat']) 
    colormap(CC); 
    %sgtitle('Koopman DMD reconstruction','FontSize',20);
    for i = 1:Length
        t = (i-1)*h;
        Z = reshape(real(DiscreteReconFun(i-1,x0)),[],Width);
        subplot(2,1,1)
        fig1 = surf(XX,YY,real(Z));
        set(fig1,'edgecolor','none');
        set(gca,'position',[0 0.5 1 0.5]);
        view(0,90);
        colorbar
        title(['Reconstructed  $x(' num2str(t) ')$'],'FontSize',20,'Interpreter','latex','Units', 'normalized', 'Position', [0.5, 0.8, 0]);
        set(gca,'YTickLabel',[],'XTickLabel',[],'XColor', 'none','YColor','none');
        subplot(2,1,2)
        ZZ = reshape(X(:,i),[],Width);
        fig1 = surf(XX,YY,real(ZZ));
        set(fig1,'edgecolor','none');
        set(gca,'position',[0 0 1 0.5]);
        view(0,90);
        colorbar
        title(['True  $x(' num2str(t) ')$'],'FontSize',20,'Interpreter','latex','Units', 'normalized', 'Position', [0.5, 0.8, 0]);
        set(gca,'YTickLabel',[],'XTickLabel',[],'XColor', 'none','YColor','none');
        drawnow;
        if SaveReconVideo
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
    end
    if SaveReconVideo
        close(v);
    end
end

%% Plots for the paper
if PrintRecon
    ReconIndex = [10 55 80 113 150];
    x = 1:1:Width;
    y = 1:1:Dimension/Width;
    x0 = X(:,1);
    [XX,YY] = meshgrid(x,y);
    load([DATAPATH '/CCcool.mat']) 
    for i=1:length(ReconIndex)
        ii = ReconIndex(i);
        t = (ii-1)*h;
        Z = reshape(real(ContinuousReconFun(t,x0)),[],Width);
        fig1 = surf(XX,YY,real(Z));
        set(fig1,'edgecolor','none');
        view(0,90);
        colormap(CC);
        set(gca,'YTickLabel',[],'XTickLabel',[],'XColor', 'none','YColor','none');
        filename = ['cyl-' K '-williams-reconstruction-at-' num2str(t) '.pdf'];
        set(gca,'position',[0 0 1 1],'YTickLabel',[],'XTickLabel',[],'XColor', 'none','YColor','none');
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [10 5]);
        set(gcf, 'PaperPosition', [0 0 10 5]);
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        saveas(gcf,filename);
    end
end