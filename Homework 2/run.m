function [mu, sigma, meanPositionError, meanMahalanobisError, meanPOfZ] = ...
	run(numSteps, usePF, pauseTime, fixSeed, doPlot, alphaBetaFactors, numParticles)
% Run a localization experiment
%
% INPUTS:
%  numSteps             - number of simulation steps
%  usePF                - Flag for using the particle filter (Default: false)
%  pauseTime            - Number of seconds to pause between two simulation steps
%                         (if <= 0, pause waiting for keypress)
%  fixSeed              - Flag for fixing the seed for the random number generator (Default: false)
%  doPlot               - Flag to turn on plotting (Default: true)
%  alphaBetaFactors     - Scaling for the noise factors of the form:
%                          [Motion noise for data, Motion noise for filter, Obs noise for data, Obs noise for filter]
%  numParticles         - Number of particles for the particle filter (Default: 100)
%
% OUTPUTS:
%  mu                   - Final state mean [x, y, theta]
%  sigma                - Final state covariance
%  meanPositionError    - Average position error
%  meanMahalanobisError - Average mahalanobis error (position error weighted by inverse co-variance)
%  meanPOfZ             - Mean of observation likelihood

% State/Observation/Control representation
%State: [x,y,theta];
%Observation: [bearing to landmark, landmark ID];
%Control: [drot1,trans,drot2];

% assume use of EKF unless usePF is present and true
if nargin < 2
  usePF = false;
end

% if pauseTime <= 0, pause waiting for keypress
if nargin < 3
  pauseTime = 0.2;
end

if nargin < 4
  fixSeed = false;
end

if nargin < 5
  doPlot = true;
end

if nargin < 6
  alphaBetaFactors = [1.0 1.0 1.0 1.0];
end

if nargin < 7
  numParticles = 100;
end

%--------------------------------------------------------------
% Graphics
%--------------------------------------------------------------
% Plot colors:
% Noise free path          - green
% Real path (ground truth) - blue
% Filter result            - red
NOISEFREE_PATH_COL = 'g';
REAL_PATH_COL = 'b';
FILTER_PATH_COL = 'r';

% Figure IDs
GLOBAL_FIGURE = 1;
ERROR_FIGURE = 2;

%--------------------------------------------------------------
% Initializations
%--------------------------------------------------------------
% Initial state mean and co-variance
initialstatemean = [180.0 50.0 0.0]';
initialsigma = [10.0, 0, 0; 0, 10.0, 0; 0, 0, 1.0];

% Motion noise (in odometry space, see p.134 in book). 
alphas = [0.05^2 0.005^2 0.1^2 0.01^2];
alphaFactorForData = alphaBetaFactors(1);
alphaFactorForFilter = alphaBetaFactors(2);

% Variance of Gaussian sensor noise
beta = deg2rad(5)^2;
betaFactorForData = alphaBetaFactors(3);
betaFactorForFilter = alphaBetaFactors(4);

% Step size between filter updates, can be less than 1.
deltaT=0.1; % check this

% consistent runs with fixed seed
if (fixSeed)
  rand('seed', 0);
  randn('seed', 0);
end

% Generate motion and sensor info consistent with noise models.
data = generateScript(initialstatemean, numSteps, alphaFactorForData*alphas, betaFactorForData*beta, deltaT);
% return values:
% data(n,:) = [realObservation', noisefreeMotion', noisefreeObservation', realRobot', noisefreeRobot'];

% Initialize ekf
mu = initialstatemean;
sigma = initialsigma;

% Initialize pf
sampleInitialParticles = true;
particles = zeros(numParticles, 3);
weights = zeros(numParticles,1);
for n = 1:numParticles
  if sampleInitialParticles
	particles(n,:) = sample(mu, sigma, 3);
  else
    particles(n,:) = mu;
  end
  weights(n) = 1.0 / numParticles; % Uniform weights
end

% print out experiment numbers:
fprintf('Alpha Factor for Data: %d \n',alphaFactorForData);
fprintf('Alpha Factor for Filter: %d \n',alphaFactorForFilter);
fprintf('Beta Factor for Data: %d \n',betaFactorForData);
fprintf('Beta Factor for Filter: %d \n',betaFactorForFilter);
if usePF
  fprintf('Using particle filter. Num particles: %d \n',numParticles);
end

%--------------------------------------------------------------
% Run the filter
%--------------------------------------------------------------

% Call ekfUpdate and pfUpdate in every iteration of this loop.
tic;
for t = 1:numSteps

  % since we don't care about the data association problem, noisefreeObs(2)
  % and realObs(2) will always be the same
  realObs = data(t,1:2)';
  markerId = int32(realObs(2));
  motionCommand = data(t,3:5)';
  noiseFreeObs = data(t,6:7)'; % not allowed for filters
  realRobot = data(t,8:10)'; % not allowed for filters, but useful for computing filter error
  noiseFreeRobot = data(t,11:13)'; % not allowed for filters
  
  % real position
  x = data(t,8);
  y = data(t,9);
  theta = data(t,10);
  
  % Plot if necessary
  if doPlot
	  figure(GLOBAL_FIGURE); clf; hold on; plotfield(markerId);
	  
	  % plot the real robot position
	  plotrobot( x, y, theta, 'k', true, 'c', 5.0);
	  
	  % draw real path and path that would result if there was no noise added to the motion
	  % this is the initial segment of these two paths (the first row of data already moves from initialstatemean)
	  plot([initialstatemean(1) data(1,8)], [initialstatemean(2) data(1,9)], 'Color', REAL_PATH_COL);
	  plot([initialstatemean(1) data(1,11)], [initialstatemean(2) data(1,12)], 'Color', NOISEFREE_PATH_COL);
	
	  % the rest of the paths
	  % draw real path
	  plot(data(1:t,8), data(1:t,9), 'Color', REAL_PATH_COL);
	  % draw noisy path
	  plot(data(1:t,11), data(1:t,12), 'Color', NOISEFREE_PATH_COL);
	
	  % indicate observed angle relative to real position
	  plot([x x+cos(theta+realObs(1))*100], [y y+sin(theta+realObs(1))*100], 'Color', REAL_PATH_COL);
  end

  %% TODO: update the filter here and output information (to the plot, for example)
  if usePF
    % [PFResults...] = pfUpdate(args)
    [particles, weights, mu, sigma, predParticles, predWeights, predMu, predSigma, pOfZ] = pfUpdate( particles, weights, numParticles, ...
						  motionCommand, alphaFactorForFilter*alphas, ...
						  realObs, betaFactorForFilter*beta, ...
						  markerId, getfieldinfo);
                      
    % You might want to plot the particles here. You could also plot co-variance matrix
    % figure(GLOBAL_FIGURE);
    % plotSamples(args)
    % plotcov2d(args)
  else
    % [EKFResults...] = ekfUpdate(args)
    [mu, sigma, predMu, predSigma, zHat, pOfZ, G, V, H, K] = ekfUpdate(mu, sigma, ...
			   motionCommand, alphaFactorForFilter*alphas, ...
			   realObs, betaFactorForFilter*beta, ...
			   markerId, getfieldinfo);
    
     % You might want to plot the mean/co-variance here
     % figure(GLOBAL_FIGURE);
     % plotcov2d(args)
  end

  %% TODO: REMOVE THIS ONCE YOU WANT TO TEST YOUR FILTERS
  mu = noiseFreeRobot;
  sigma = initialsigma;
  pOfZ = 1.0;
  disp('Currently ignoring filter values');
  %% END REMOVE

  % Concatenate results
  muList(t,:) = mu;
  pOfZList(t) = pOfZ;
  
  % Plot filter path, current mean and co-variance
  if (doPlot)
    figure(GLOBAL_FIGURE);
    plot(muList(1:t,1), muList(1:t,2), 'Color', FILTER_PATH_COL);
  end

  %% TODO?: you may want to plot and evaluate additional filter results here
  % example of a second plot (feel free to remove):
  % positionDriftExample = sqrt(sumsq( (data(1:t,8:9) - data(1:t,11:12)), 2));
  % figure(ERROR_FIGURE); clf; plot(positionDriftExample);

  errors(t,:) = mu - realRobot;
  errors(t,3) = minimizedAngle(errors(t,3));
  positionErrors(t) = norm(errors(t,1:2));
  condNumber = cond(sigma);
  if condNumber > 1e12
    disp('badly conditioned sigma (setting to identity):');
    sigma
    condNumber
    sigma = eye(3);
    fflush(stdout);
  end
  mahalanobis(t) = errors(t,:)*inv(sigma)*errors(t,:)';

  if pauseTime > 0
    pause(pauseTime);
  else
    disp('paused...');
    pause();
  end
end
toc;

% Display results
meanPositionError = mean(positionErrors)
meanMahalanobisError = mean(mahalanobis)
ANEES = meanMahalanobisError / 3
meanPOfZ = mean(pOfZList)
disp('');

% Save figure as png file
figure(GLOBAL_FIGURE); print('filter-path.png', '-dpng');

