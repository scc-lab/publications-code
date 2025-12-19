function [resT,msh,longestPeriod] = cantileverBeamUndampedTransient(F,zeta)
%% Dynamics of Damped Cantilever Beam
% This file contains a slighty modified version of Code taken from 
% https://www.mathworks.com/help/pde/ug/dynamics-of-a-damped-cantilever-beam.html
% 
% This example shows how to include damping in the transient analysis of a simple 
% cantilever beam.
% 
% The damping model is basic viscous damping distributed uniformly through the 
% volume of the beam. The beam is deformed by applying an external load at the 
% tip of the beam and then released at time $t = 0$. This example does not use 
% any additional loading, so the displacement of the beam decreases as a function 
% of time due to the damping. The example uses plane-stress modal, static, and 
% transient analysis models in its three-step workflow:
%% 
% # Perform modal analysis to compute the fundamental frequency of the beam 
% and to speed up computations for the transient analysis.
% # Find the static solution of the beam with a vertical load at the tip to 
% use as an initial condition for a transient model.
% # Perform the transient analysis with and without damping.
%% 
% Damping is typically expressed as a percentage of critical damping, $\xi$, 
% for a selected vibration frequency. This example uses $\xi = 0.03$, which is 
% three percent of critical damping.
% 
% The example specifies values of parameters using the imperial system of units. 
% You can replace them with values specified in the metric system. If you do so, 
% ensure that you specify all values throughout the example using the same system.
%% Modal Analysis
% Create a modal analysis model for a plane-stress problem.

modelM = createpde('structural','modal-planestress');
%% 
% Create the geometry and include it in the model. Suppose, the beam is 5 inches 
% long and 0.1 inches thick.

width = 5; 
height = 0.1;

gdm = [3;4;0;width;width;0;0;0;height;height];
g = decsg(gdm,'S1',('S1')');
geometryFromEdges(modelM,g);
%% 
% Define a maximum element size so that there are five elements through the 
% beam thickness. Generate a mesh.

hmax = height/5;
msh = generateMesh(modelM,'Hmax',hmax);
%% 
% Specify the Young's modulus, Poisson's ratio, and mass density of steel.

E = 3.0e7; 
nu = 0.3; 
rho = 0.3/386;
structuralProperties(modelM,'YoungsModulus',E, ...
                            'PoissonsRatio',nu, ...
                            'MassDensity',rho);
%% 
% Specify that the left edge of the beam is a fixed boundary.

structuralBC(modelM,'Edge',4,'Constraint','fixed');
%% 
% Solve the problem for the frequency range from |0| to |1e5|. The recommended 
% approach is to use a value that is slightly smaller than the expected lowest 
% frequency. Thus, use |-0.1| instead of |0|.

res = solve(modelM,'FrequencyRange',[-0.1,1e5]');
%% 
% By default, the solver returns circular frequencies.

modeID = 1:numel(res.NaturalFrequencies);
%% 
% Express the resulting frequencies in Hz by dividing them by |2π|_._ Display 
% the frequencies in a table.

tmodalResults = table(modeID.',res.NaturalFrequencies/(2*pi));
tmodalResults.Properties.VariableNames = {'Mode','Frequency'};
%disp(tmodalResults)
%% 
% Compute the analytical fundamental frequency (Hz) using the beam theory.

I = height^3/12;
freqAnalytical = 3.516*sqrt(E*I/(width^4*rho*height))/(2*pi);
%% 
% Compare the analytical result with the numerical result.

freqNumerical = res.NaturalFrequencies(1)/(2*pi);
%% 
% Compute the period corresponding to the lowest vibration mode.

longestPeriod = 1/freqNumerical;
%% Initial Displacement from Static Solution
% The beam is deformed by applying an external load at its tip and then released 
% at time $t = 0$. Find the initial condition for the transient analysis by using 
% the static solution of the beam with a vertical load at the tip.
% 
% Create a static plane-stress model.

modelS = createpde('structural','static-planestress');
%% 
% Use the same geometry and mesh that you used for the modal analysis.

geometryFromEdges(modelS,g);
modelS.Mesh = msh;
%% 
% Specify the same values for the Young's modulus, Poisson's ratio, and mass 
% density of the material.

structuralProperties(modelS,'YoungsModulus',E, ...
                           'PoissonsRatio',nu, ...
                           'MassDensity',rho);
                       
%% 
% Specify the same constraint on the left end of the beam.

structuralBC(modelS,'Edge',4,'Constraint','fixed');
%% 
% Apply the static vertical load on the right side of the beam.

structuralBoundaryLoad(modelS,'Edge',2,'SurfaceTraction',[0;F]);
%% 
% Solve the static model. The resulting static solution serves as an initial 
% condition for transient analysis.

Rstatic = solve(modelS);
%% Transient Analysis
% Perform the transient analysis of the cantilever beam with and without damping. 
% Use the modal superposition method to speed up computations.
% 
% Create a transient plane-stress model.

modelT = createpde('structural','transient-planestress');
%% 
% Use the same geometry and mesh that you used for the modal analysis.

geometryFromEdges(modelT,g);
modelT.Mesh = msh;
%% 
% Specify the same values for the Young's modulus, Poisson's ratio, and mass 
% density of the material.

structuralProperties(modelT,'YoungsModulus',E, ...
                            'PoissonsRatio',nu, ...
                            'MassDensity',rho);
                       
%% 
% Specify the same constraint on the left end of the beam.

structuralBC(modelT,'Edge',4,'Constraint','fixed');
%% 
% Specify the initial condition by using the static solution.

structuralIC(modelT,Rstatic);
%% 
% Solve the undamped transient model for three full periods corresponding to 
% the lowest vibration mode.

tlist = 0:longestPeriod/100:3*longestPeriod;
resT = solve(modelT,tlist,'ModalResults',res);
%% 
% Interpolate the displacement at the tip of the beam.

intrpUt = interpolateDisplacement(resT,[5;0.05]);
%% 
% Now solve the model with no damping.

%zeta = 0.0;
omega = 2*pi*freqNumerical;
structuralDamping(modelT,'Zeta',zeta);
resT = solve(modelT,tlist,'ModalResults',res);
