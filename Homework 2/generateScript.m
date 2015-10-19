%generateScript: simulates the trajectory of the robot using square
%                path given by generate motion
%
%data=generateScript(initialstatemean,numSteps,alphas,betas)
%     generates data of the form 
%     [realObservation', noisefreeMotion', noisefreeObservation',
%     realRobot', noisefreeRobot']
%
%realObservation and noisefreeMotion is the only data available to
%filter.  All other data for debugging/display purposes.
%
%alphas are the 4-d noise for robot motion
%beta: noise for observations
%

function data = generateScript(initialstatemean, numSteps, alphas, ...
			       beta, deltaT)

%--------------------------------------------------------------
% Initializations
%--------------------------------------------------------------

global FIELDINFO;
FIELDINFO = getfieldinfo; % retreive information about the soccer field

realRobot = initialstatemean;
noisefreeRobot = initialstatemean;

for n = 1:numSteps
  % --------------------------------------------
  % Simulate motion
  % --------------------------------------------

  t=n*deltaT;
  noisefreeMotion = generateMotion(t,deltaT);
  
  % Shift real robot
%   prevNoisefreeRobot = noisefreeRobot;
  noisefreeRobot = prediction( noisefreeRobot, noisefreeMotion);
  
  % Move robot
  realRobot = sampleOdometry( noisefreeMotion, realRobot, alphas);
  
  %--------------------------------------------------------------
  % Simulate observation
  %--------------------------------------------------------------

  % n / 2 causes each landmark to be viewed twice
  markerId = mod( floor( n / 2), FIELDINFO.NUM_MARKERS) + 1;
  
  noisefreeObservation = observation( realRobot, FIELDINFO, markerId);

  % Observation noise
  observationNoise = sample( 0.0, beta, 1);
  realObservation = noisefreeObservation;
  realObservation(1) = realObservation(1) + observationNoise;
  
  data(n,:) = [realObservation', noisefreeMotion', noisefreeObservation', realRobot', noisefreeRobot'];
end
