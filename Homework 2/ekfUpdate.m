function [mu, sigma, predMu, predSigma, zHat, pOfZ, G, V, H, K] = ...
  	ekfUpdate( mu, sigma, ...
			   u, filterAlphas, ...
			   z, filterBeta, ...
			   markerId, FIELDINFO)
% EKF Update step
%
% INPUT:
%  mu           - state mean at time t-1
%  sigma        - state covariance at time t-1
%  u            - control at time t
%  filterAlphas - motion model noise
%  z            - observation at time t
%  filterBeta   - observation noise variance
%  markerId     - ID of the marker we see at time t (get observation from)
%  FIELDINFO    - Field details for observation function
%
% OUTPUTS:
%  mu           - state mean at time t
%  sigma        - state covar at time t 
%  predMu       - state mean from EKF prediction step
%  predSigma    - state covar from EKF prediction step
%  zHat         - expected observation
%  pOfZ         - observation likelihood
%  G            - Jacobian of dynamics w.r.t state_t-1
%  V            - Jacobian of dynamics w.r.t control_t
%  H            - Observation jacobian
%  K            - Kalman gain

% TODO: Remove this line
[predMu, predSigma, zHat, pOfZ, G, V, H, K] = deal([]);

% You should almost certainly leave the header the way it is.

% Init vars
Delta_rot1 = u(1);
Delta_trans = u(2);
Delta_rot2 = u(3);
Pos_prev_x = mu(1);
Pos_prev_y = mu(2);
Pos_prev_theta = mu(3);

% --------------------------------------------
% Prediction step
% --------------------------------------------

% some stuff here


%--------------------------------------------------------------
% Correction step
%--------------------------------------------------------------

% Compute expected observation and Jacobian

% Innovation / residual covariance

% Residual

% Likelihood

% Kalman gain

% Correction





   
