%% Paper title "Safe Model-Based Reinforcement Learning for Systems with Parametric Uncertainties" 
%% arXiv:2007.12666
%% Coded by S M Nahid Mahmud, MS Grad Student, Oklahoma State University.
%% nahid.mahmud@okstate.edu
%% Four-state dynamical system (Section 7.2 of the journal paper)

clear all; 
clc
 
%-------------------------------------------------------------------------%
%----------------- Provide All Bounds for Problem ------------------------%
%-------------------------------------------------------------------------%
a1 = -7; A1 = 5; %Barrier Boundaries for x1, x1-coordinate values will be remained within (a1,A1) 
a2 = -7; A2 = 5; %Barrier Boundaries for x2, x2-coordinate values will be remained within (a2,A2) 
a3 = -5; A3 = 7; %Barrier Boundaries for x3, x3-coordinate values will be remained within (a3,A3) 
a4 = -5; A4 = 7; %Barrier Boundaries for x4, x4-coordinate values will be remained within (a4,A4) 


% Initial x-coordinate values, it has to be within (a,A).    
 x1_now0 = -5;
 x2_now0 = -5;
 x3_now0 = 5;
 x4_now0 = 5;
    
%transforming initial x-coordinates to initial s-coordinates. 
s1_now = (log ((A1/a1)*((a1-x1_now0))/(A1-x1_now0)));
s2_now = (log ((A2/a2)*((a2-x2_now0))/(A2-x2_now0))); 
s3_now = (log ((A3/a3)*((a3-x3_now0))/(A3-x3_now0)));
s4_now = (log ((A4/a4)*((a4-x4_now0))/(A4-x4_now0))); 

% Giving upper boundary to x coordinates for this simulation purpose
x1_max = A1-0.0001;
x2_max = A2-0.0001;  
x3_max = A3-0.0001;
x4_max = A4-0.0001; 
    
% Giving lower boundary to x coordinates for this simulation purpose    
x1_min = a1+0.0001;
x2_min = a2+0.0001;
x3_min = a3+0.0001;
x4_min = a4+0.0001;
    
% Giving upper boundary to s coordinates for this simulation purpose 
s1_max = (log ((A1/a1)*((a1-x1_max))/(A1-x1_max)));
s2_max = (log ((A2/a2)*((a2-x2_max))/(A2-x2_max))); 
s3_max = (log ((A3/a3)*((a3-x3_max))/(A3-x3_max)));
s4_max = (log ((A4/a4)*((a4-x4_max))/(A4-x4_max))); 
  
% Giving lower boundary to s coordinates for this simulation purpose      
s1_min = (log ((A1/a1)*((a1-x1_min))/(A1-x1_min)));
s2_min = (log ((A2/a2)*((a2-x2_min))/(A2-x2_min))); 
s3_min = (log ((A3/a3)*((a3-x3_min))/(A3-x3_min)));
s4_min = (log ((A4/a4)*((a4-x4_min))/(A4-x4_min))); 


t0 = 0; %Lower time bound
tf = 200; %Upper time bound

%Starting s-coordinate position
x0 = [s1_now s2_now s3_now s4_now];
%End s-coordinate position
xf = [0 0 0 0];

%Upper bound of ending position
xfmax = inf*[1,1,1,1];
%Lower bound of ending position
xfmin = -inf*[1,1,1,1];

%Upper bound of starting position
xMax = [s1_max s2_max s3_max s4_max ];
%Lower bound of starting position
xMin = [s1_min s2_min s3_min s4_min];




%Lower bound of controller
uMin = [-10 -10];
%Upper bound of controller
uMax = [10 10];

 
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
guess.phase.control  = [0 0;0 0] 
guess.phase.integral = 0;
 
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp-LiuRao-Legendre';
mesh.tolerance       = 1e-7;
mesh.maxiterations   = 4500;
mesh.colpointsmin    = 50;
mesh.colpointsmax    = 200;
mesh.phase.colpoints = 4*ones(1,10);
mesh.phase.fraction  = 0.1*ones(1,10);
 
%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name                           = 'four-state';
setup.functions.continuous           = @continuous;
setup.functions.endpoint             = @endpoint;
setup.displaylevel                   = 2;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.nlp.solver                     = 'ipopt';
setup.nlp.snoptoptions.tolerance     = 1e-8;
setup.nlp.snoptoptions.maxiterations = 20000;
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.tolerance     = 1e-8;
setup.derivatives.supplier           = 'SparseCD';
setup.derivatives.derivativelevel    = 'second';
setup.method                         = 'RPM-Differentiation';
 
%-------------------------------------------------------------------------%
%---------------------- Solve Problem Using GPOPS2 -----------------------%
%-------------------------------------------------------------------------%
tic
output = gpops2(setup);
toc

