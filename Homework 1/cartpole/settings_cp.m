%Set random seed
rng('default');

% include some paths
try
  rd = 'cartpole/';
  addpath([rd 'base'], [rd 'util'], [rd 'control']);
catch
end

% 1. Define state and important indices

% 1a. Full state representation (including all augmentations)
%
%  1  x          cart position
%  2  v          cart velocity
%  3  dtheta     angular velocity
%  4  theta      angle of the pendulum
%  5  sin(theta) complex representation ...
%  6  cos(theta) of theta
%  7  u          force applied to cart
%

% 1b. Important indices
% odei  indicies for the ode solver
% dyno  indicies for the output from the dynamics model and indicies to loss
% angi  indicies for variables treated as angles (using sin/cos representation)
% poli  indicies for the inputs to the policy

% 2. Set up the scenario
dt = 0.1;                         % [s] sampling time
totT = 3.0;                           % [s] initial prediction horizon time
H = ceil(totT/dt);                    % prediction steps (optimization horizon)
mu0 = [0 0 0 0]';                  % initial state mean
S0 = diag([0.1 0.01 0.01 pi/2].^2);   % initial state covariance

% 3. Some setup for the full dynamics simulator
plant.dynamics = @dynamics_cp;
plant.noise = diag(ones(1,4)*0.01.^2);            % measurement noise
plant.dt = dt;
plant.ctrl = @zoh;                                % controler is zero order hold
plant.odei = [1 2 3 4];            % varibles for the ode solver
plant.angi = [4];                  % angle variables
plant.poli = [1 2 3 5 6];          % variables that serve as inputs to the policy
plant.dyno = [1 2 3 4];            % variables to be predicted (and known to loss)

% 4. Pre-set policy
dataset = load('cartPole_10_H30.mat');
dataset.policy.fcn = @(policy,m,s)conCat(@conpred,@gSat,policy,m,s);
policy.maxU = 10;                                         % max. amplitude of 

% % Old hyper-parameters
% kernelLengthScale = [241.5374  253.5284  267.9354  258.0687
%                      28.1628  215.8184  213.2465  220.1273
%                      38.9766    7.6765    7.4598   14.1591
%                       2.0037    1.4048    0.9705    1.6362
%                      29.1488    0.8736    0.7849    1.6199
%                     195.2348   27.3041   24.2517   65.3293];
% kernelScaleFactor = [1.0236    1.5798    3.8723    0.8942];
% noiseSigma        = [0.0157    0.0155    0.0411    0.0153];
% New ones tuned for:
% Tuned for random policy + learned policy (every 4 steps), 10 control, 
% 0.1 dt, 3 secs,  mu = 0, S varies based on policy

% 5. Setup Figures for plotting
fg = figure(101); hold on;
set(fg, 'Position', [50, 100, 1200, 550]);
pltHdls = zeros(1,6);
pltHdls(1) = subplot('Position', [0.05,0.5,0.4,0.5]);
pltHdls(2) = subplot('Position', [0.55,0.5,0.4,0.5]);
pltHdls(3) = subplot('Position', [0.095,0.33,0.33,0.15]);
pltHdls(4) = subplot('Position', [0.595,0.33,0.33,0.15]);
pltHdls(5) = subplot('Position', [0.095,0.1,0.33,0.15]);
pltHdls(6) = subplot('Position', [0.595,0.1,0.33,0.15]);

% Title and legend text for bottom figures
TT = subplot('Position',[0.5,0.55,0.01,0.01]);
text(-14,0,'Comparison of GP Predictions (\Delta{State}) vs Full Dynamics', ...
  'fontweight','bold','FontName','Times New Roman','fontsize',14);
axis off;

TT1 = subplot('Position',[0.49,0.31,0.05,0.04]); 
text(0,0,'Full-Dynamics','fontweight','bold','FontName','Times New Roman','fontsize',12);
text(0,-1,'GP-Dynamics','fontweight','bold','FontName','Times New Roman','fontsize',12);
axis off; 

TT2 = subplot('Position',[0.43,0.265,0.05,0.04]); hold on;
plot([0 -0.8],[0, 0],'g-','linewidth',3)
plot([0 -0.8],[-1, -1],'r-','linewidth',3)
axis off; hold off;

predTxt = {'Cart Vel (m/s)', 'Cart Acc (m/s^2)', 'Ang Acc (rad/s^2)', 'Ang Vel (rad/s)'};
