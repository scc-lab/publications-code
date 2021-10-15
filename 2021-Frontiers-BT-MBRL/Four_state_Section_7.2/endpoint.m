%% Paper title "Safe Model-Based Reinforcement Learning for Systems with Parametric Uncertainties" 
%% arXiv:2007.12666
%% Coded by S M Nahid Mahmud, MS Grad Student, Oklahoma State University.
%% nahid.mahmud@okstate.edu
%% Four-state dynamical system (Section 7.2 of the journal paper)

function output = endpoint(input)
q  = input.phase.integral;
output.objective = q;












