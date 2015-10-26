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

% Init vars
Delta_rot1 = u(1);
Delta_trans = u(2);
Delta_rot2 = u(3);
Pos_prev_x = mu(1);
Pos_prev_y = mu(2);
Pos_prev_theta = mu(3);
g=@prediction;

% --------------------------------------------
% Prediction step
% --------------------------------------------

% some stuff here
predMu  = g(mu, u);

% finite difference
h = 1e-4;

%% Produce G
muDelx = [(Pos_prev_x + h)
            Pos_prev_y
            Pos_prev_theta];
muDely = [(Pos_prev_x)
    Pos_prev_y + h
    Pos_prev_theta];
muDeltheta = [Pos_prev_x
    Pos_prev_y
    Pos_prev_theta + h];
predMuDelX = prediction(muDelx, u);
predMuDelY = prediction(muDely, u);
predMuDelTheta = prediction(muDeltheta, u);

% example for first col, changing x
x_0 = predMu(1);
x_1 = predMuDelX(1);

y_0 = predMu(2);
y_1 = predMuDelX(2);

theta_0 = predMu(3);
theta_1 = predMuDelX(3);

dx = x_1 - x_0;
dy = y_1 - y_0;
dtheta = theta_1 - theta_0;

G1 = [dx/h
      dy/h
      dtheta/h];
  
  
% example for second col, changing y
x_0 = predMu(1);
x_1 = predMuDelY(1);

y_0 = predMu(2);
y_1 = predMuDelY(2);

theta_0 = predMu(3);
theta_1 = predMuDelY(3);

dx = x_1 - x_0;
dy = y_1 - y_0;
dtheta = theta_1 - theta_0;

G2 = [dx/h
      dy/h
      dtheta/h];
  
% example for third col, changing theta
x_0 = predMu(1);
x_1 = predMuDelTheta(1);

y_0 = predMu(2);
y_1 = predMuDelTheta(2);

theta_0 = predMu(3);
theta_1 = predMuDelTheta(3);

dx = x_1 - x_0;
dy = y_1 - y_0;
dtheta = theta_1 - theta_0;

G3 = [dx/h
      dy/h
      dtheta/h];

G = horzcat(G1, G2);
G = horzcat(G, G3);


%% Produce V

uDelx = [Delta_rot1 + h
         Delta_trans
         Delta_rot2];
uDely = [Delta_rot1
         Delta_trans + h
         Delta_rot2];
uDeltheta = [Delta_rot1
         Delta_trans
         Delta_rot2 + h];
preduDelX = prediction(mu, uDelx);
preduDelY = prediction(mu, uDely);
preduDelTheta = prediction(mu, uDeltheta);

% example for first col, changing x
x_0 = predMu(1);
x_1 = preduDelX(1);

y_0 = predMu(2);
y_1 = preduDelX(2);

theta_0 = predMu(3);
theta_1 = preduDelX(3);

dx = x_1 - x_0;
dy = y_1 - y_0;
dtheta = theta_1 - theta_0;

V1 = [dx/h
      dy/h
      dtheta/h];

% second col, changing y control
x_0 = predMu(1);
x_1 = preduDelY(1);

y_0 = predMu(2);
y_1 = preduDelY(2);

theta_0 = predMu(3);
theta_1 = preduDelY(3);

dx = x_1 - x_0;
dy = y_1 - y_0;
dtheta = theta_1 - theta_0;

V2 = [dx/h
      dy/h
      dtheta/h];
  
  
% third col, changing theta control
x_0 = predMu(1);
x_1 = preduDelTheta(1);

y_0 = predMu(2);
y_1 = preduDelTheta(2);

theta_0 = predMu(3);
theta_1 = preduDelTheta(3);

dx = x_1 - x_0;
dy = y_1 - y_0;
dtheta = theta_1 - theta_0;

V3 = [dx/h
      dy/h
      dtheta/h];

% Combine to form V
V = horzcat(V1, V2);
V = horzcat(V, V3);


%% Produce SigmaHat
M_t = diag(noiseFromMotion(u, filterAlphas));

sigmaHat = G * sigma * transpose(G) + V * M_t * transpose(V);
predSigma = sigmaHat;
%% --------------------------------------------------------------
% Correction step
%--------------------------------------------------------------

%% Compute expected observation and Jacobian

zHat = observation(mu, FIELDINFO, markerId);

zHatDelx = observation(muDelx, FIELDINFO, markerId);
zHatDely = observation(muDely, FIELDINFO, markerId);
zHatDeltheta = observation(muDeltheta, FIELDINFO, markerId);

b_0 = zHat(1);
b_1_Delx = zHatDelx(1);
b_1_Dely = zHatDely(1);
b_1_Deltheta = zHatDeltheta(1);


delbDelx = (b_1_Delx - b_0) / h;
delbDely = (b_1_Dely - b_0) / h;
delbDeltheta = (b_1_Deltheta - b_0) / h;

H = [ delbDelx delbDely delbDeltheta ];

Q =  [  [filterBeta(1).^ 2 0
            0               0]  ];
Q = filterBeta(1)^2;

S = H * sigmaHat * transpose(H) + Q;

   
%% Innovation / residual covariance

innovation = z(1) - zHat(1);

%% Residual

%% Likelihood

%% Kalman gain

K = sigmaHat * transpose(H) * inv(S);


%% Correction

mu = predMu + K * innovation; % TODO: add Kalman gain * innovation
I = eye(3);
sigma = (I - K * H) * sigmaHat; % TODO: this is not done
%sigma = sigmaHat;
%mu = predMu;

pOfZ = sqrt(det(2 * pi * S)) * exp(-0.5 * transpose(innovation) * inv(S) * innovation);
testLine = 5;


   
