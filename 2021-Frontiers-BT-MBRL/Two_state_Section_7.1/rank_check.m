%% Paper title "Safe Model-Based Reinforcement Learning for Systems with Parametric Uncertainties" 
%% arXiv:2007.12666
%% Coded by S M Nahid Mahmud, MS Grad Student, Oklahoma State University.
%% nahid.mahmud@okstate.edu
%% Two-state dynamical system (Section 7.1 of the journal paper)

%% This file is to check for the time require for the ScriptYf matrix to fulfill the full rank condition.  
%% Row 15 to 30 of z means ScriptYf matrix elements
Rank = load('States_RL_BF_FCL_manipulator.mat');
A = [Rank.z(:,15),Rank.z(:,19),Rank.z(:,23),Rank.z(:,27)];
B = [Rank.z(:,16),Rank.z(:,20),Rank.z(:,24),Rank.z(:,28)];
C = [Rank.z(:,17),Rank.z(:,21),Rank.z(:,25),Rank.z(:,29)];
D = [Rank.z(:,18),Rank.z(:,22),Rank.z(:,26),Rank.z(:,30)];
for i = 1: 1797965 %this 1797965 number is from the total iteration number we got from the simulation.
%This 1797965 should be changed with every simulation run
E = [A(i,:);B(i,:);C(i,:);D(i,:)]
F = rank(E) 
if F == 4
    i
       disp(stop)
end 
end 