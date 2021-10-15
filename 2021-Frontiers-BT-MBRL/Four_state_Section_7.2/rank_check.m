%% Paper title "Safe Model-Based Reinforcement Learning for Systems with Parametric Uncertainties" 
%% arXiv:2007.12666
%% Coded by S M Nahid Mahmud, MS Grad Student, Oklahoma State University.
%% nahid.mahmud@okstate.edu
%% Four-state dynamical system (Section 7.2 of the journal paper)
%% This file is to check for the time require for the ScriptYf matrix to fulfill the full rank condition.  
%% Row 25 to 40 of z means ScriptYf matrix elements

Rank = load('States_RL_BF_FCL_manipulator.mat');
A = [Rank.z(:,25),Rank.z(:,29),Rank.z(:,33),Rank.z(:,37)];
B = [Rank.z(:,26),Rank.z(:,30),Rank.z(:,34),Rank.z(:,38)];
C = [Rank.z(:,27),Rank.z(:,31),Rank.z(:,35),Rank.z(:,39)];
D = [Rank.z(:,28),Rank.z(:,32),Rank.z(:,36),Rank.z(:,40)];
for i = 1: 1797965 %this 1797965 number is from the total iteration number we got from the simulation.
%This 1797965 should be changed with every simulation run
E = [A(i,:);B(i,:);C(i,:);D(i,:)]
F = rank(E) 
if F == 4
    i
       disp(stop)
end 


end 
