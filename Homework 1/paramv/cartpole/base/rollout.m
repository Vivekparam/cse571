%% rollout.m
% *Summary:* Generate a state trajectory using an ODE solver (and any additional 
% dynamics) from a particular initial state by applying either a particular 
% policy or random actions.
%
%   function [x y] = rollout(start, policy, H, plant)
%
% *Input arguments:*
%   
%   start       vector containing initial states (without controls)   [nX  x  1]
%   policy      policy structure
%     .fcn        policy function
%     .p          parameter structure (if empty: use random actions)
%     .maxU       vector of control input saturation values           [nU  x  1]
%   H           rollout horizon in steps
%   plant       the dynamical system structure
%     .poli       indices for states passed to the policy
%     .dyno       indices for states passed to cost
%     .odei       indices for states passed to the ode solver
%
% *Output arguments:*
%
%   x          matrix of observed states                           [H   x nX+nU]
%   y          matrix of corresponding observed successor states   [H   x   nX ]                          [H+1 x   nX ]
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
%
% Modified by Arunkumar Byravan for CSE 571: Probabilistic Robotics 
% (Autumn 2015) at  University of Washington
%
%% High-Level Steps
%
% # Compute control signal $u$ from state $x$:
% either apply policy or random actions
% # Simulate the true dynamics for one time step using the current pair $(x,u)$
% # Apply random noise to the successor state
% # Repeat until end of horizon

function [x, y] = rollout(start, policy, H, plant)
%% Code

odei = plant.odei; poli = plant.poli; dyno = plant.dyno; angi = plant.angi;
simi = sort(odei);
nX = length(simi); nU = length(policy.maxU); nA = length(angi);

state(simi) = start;      % initializations
x = zeros(H+1, nX+2*nA);
x(1,simi) = start' + randn(size(simi))*chol(plant.noise);
u = zeros(H, nU);
next = zeros(1,length(simi));

for i = 1:H % --------------------------------------------- generate trajectory
  s = x(i,dyno)'; sa = gTrig(s, zeros(length(s)), angi); s = [s; sa];
  x(i,end-2*nA+1:end) = s(end-2*nA+1:end);
  
  % 1. Apply policy ... or random actions --------------------------------------
  if isfield(policy, 'fcn')
    u(i,:) = policy.fcn(policy,s(poli),zeros(length(poli)));
  else
    u(i,:) = policy.maxU.*(2*rand(1,nU)-1);
  end

  % 2. Simulate dynamics -------------------------------------------------------
  next(odei) = simulate(state(odei), u(i,:), plant);

  % 3. Augment state and randomize ---------------------------------------------
  state(simi) = next(simi);
  x(i+1,simi) = state(simi) + randn(size(simi))*chol(plant.noise);
end

y = x(2:H+1,1:nX); x = [x(1:H,:) u(1:H,:)]; 