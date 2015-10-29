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
%  markerId      - ID of the marker we see at time t (get observation from)
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
% [predMu, predSigma, zHat, pOfZ, G, V, H, K] = deal([]);

% You should almost certainly leave the header the way it is.

% --------------------------------------------
% Prediction step
% --------------------------------------------

% some stuff here
predMu  = prediction(mu, u);

% finite difference
h = 1e-7;


%% Produce G
G = zeros(3, 3);
for i = 1:3
    temp = mu;
    temp(i) = temp(i) + h;
    predMuDel = prediction(temp, u);
    G(:, i) = predMuDel - predMu;
end

G = G / h;

%% Produce V
V = zeros(3, 3);

for i = 1:3
    temp = u;
    temp(i) = temp(i) + h;
    predMuDel = prediction(mu, temp);
    V(:, i) = predMuDel - predMu;
end

V = V / h;

%% Produce SigmaHat
M_t = diag(noiseFromMotion(u, filterAlphas));

sigmaHat = G * sigma * transpose(G) + V * M_t * transpose(V);
predSigma = sigmaHat;
%% --------------------------------------------------------------
% Correction step
%--------------------------------------------------------------

%% Compute expected observation and Jacobian

zHat = observation(predMu, FIELDINFO, markerId);

H = zeros(1, 3);

for i = 1:3
   temp = predMu;
   temp(i) = temp(i) + h;
   predZHatDel = observation(temp, FIELDINFO, markerId);
   H(i) = minimizedAngle(predZHatDel(1) - zHat(1));
end

H = H / h;

Q = filterBeta(1);

S = H * sigmaHat * transpose(H) + Q;

   
%% Innovation / residual covariance

innovation = minimizedAngle(z(1) - zHat(1));

%% Residual

%% Likelihood

%% Kalman gain

K = sigmaHat * transpose(H) * inv(S);


%% Correction

mu = predMu + K * innovation;
I = eye(3);
sigma = (I - K * H) * sigmaHat;

pOfZ = likelihood(S, innovation); %sqrt(det(2 * pi * S)) * exp(-0.5 * transpose(innovation) * inv(S) * innovation);



   
