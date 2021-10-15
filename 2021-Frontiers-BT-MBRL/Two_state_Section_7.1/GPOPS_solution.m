%% Paper title "Safe Model-Based Reinforcement Learning for Systems with Parametric Uncertainties" 
%% arXiv:2007.12666
%% Coded by S M Nahid Mahmud, MS Grad Student, Oklahoma State University.
%% nahid.mahmud@okstate.edu
%% Two-state dynamical system (Section 7.1 of the journal paper)

clear all; 
clc
 
%-------------------------------------------------------------------------%
%----------------- Provide All Bounds for Problem ------------------------%
%-------------------------------------------------------------------------%

a1 = -7; A1 = 5; %Barrier Boundaries for x1, x1-coordinate values will be remained within (a1,A1) 
a2 = -5; A2 = 7; %Barrier Boundaries for x2, x2-coordinate values will be remained within (a2,A2) 


% Initial x-coordinate values, it has to be within (a,A). 
  
x1_now0 = -6.5;
x2_now0 = 6.5;
 
%transforming initial x-coordinates to initial s-coordinates.  
    
s1_now = (log ((A1/a1)*((a1-x1_now0))/(A1-x1_now0)));
s2_now = (log ((A2/a2)*((a2-x2_now0))/(A2-x2_now0))); 

% Giving upper boundary to x coordinates for this simulation purpose
x1_max = A1-0.0001;
x2_max = A2-0.0001;  

    
% Giving lower boundary to x coordinates for this simulation purpose    
x1_min = a1+0.0001;
x2_min = a2+0.0001;

    
% Giving upper boundary to s coordinates for this simulation purpose 
s1_max = (log ((A1/a1)*((a1-x1_max))/(A1-x1_max)));
s2_max = (log ((A2/a2)*((a2-x2_max))/(A2-x2_max))); 

  
% Giving lower boundary to s coordinates for this simulation purpose      
s1_min = (log ((A1/a1)*((a1-x1_min))/(A1-x1_min)));
s2_min = (log ((A2/a2)*((a2-x2_min))/(A2-x2_min))); 

t0 = 0; %Lower time bound
tf = 20; %Upper time bound


%Starting s-coordinate position
x0 = [s1_now s2_now];
%End s-coordinate position
xf = [0 0];

%Upper bound of ending position
xfmax = [0,0];
%Lower bound of ending position
xfmin = [0,0];

%Upper bound of starting position
xMax = [s1_max s2_max];
%Lower bound of starting position
xMin = [s1_min s2_min];

%Lower bound of controller
uMin = -100;
%Upper bound of controller
uMax = +100;
 
%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf;
bounds.phase.finaltime.upper = tf;
bounds.phase.initialstate.lower = x0; 
bounds.phase.initialstate.upper = x0;
bounds.phase.finalstate.lower = xfmin; 
bounds.phase.finalstate.upper = xfmax;
bounds.phase.control.lower = uMin; 
bounds.phase.control.upper = uMax;
bounds.phase.state.lower = xMin; 
bounds.phase.state.upper = xMax;
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 100000;
 
%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time     = [t0; tf]; 
guess.phase.state    = [x0; xf];
guess.phase.control  = [0; 0];
guess.phase.integral = 0;
 
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp-LiuRao-Legendre';
mesh.tolerance       = 1e-7;
mesh.maxiterations   = 4500;
mesh.colpointsmin    = 200;
mesh.colpointsmax    = 1400;
mesh.phase.colpoints = 4*ones(1,10);
mesh.phase.fraction  = 0.1*ones(1,10);
 
%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name                           = 'two-state-system';
setup.functions.continuous           = @continuous;
setup.functions.endpoint             = @endpoint;
setup.displaylevel                   = 2;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.nlp.solver                     = 'ipopt';
setup.nlp.snoptoptions.tolerance     = 1e-12;
setup.nlp.snoptoptions.maxiterations = 20000;
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.tolerance     = 1e-10;
setup.derivatives.supplier           = 'SparseCD';
setup.derivatives.derivativelevel    = 'second';
setup.method                         = 'RPM-Differentiation';
 
%-------------------------------------------------------------------------%
%---------------------- Solve Problem Using GPOPS2 -----------------------%
%-------------------------------------------------------------------------%
tic
output = gpops2(setup);
toc

