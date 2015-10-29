function [particles, weights, mu, sigma, ...
	  predParticles, predWeights, predMu, predSigma, pOfZ] = pfUpdate( particles, weights, numParticles, ...
						  u, alphas, ...
						  z, beta, ...
						  markerId, FIELDINFO)
% PF Update step
%
% INPUTS:
%  particles    - particles at time t-1
%  sigma        - particle weights at time t-1
%  numParticles - number of particles
%  u            - control at time t
%  alphas       - motion model noise
%  z            - observation at time t
%  beta         - observation noise variance
%  markerId     - ID of the marker we see at time t (get observation from)
%  FIELDINFO    - Field details for observation function
%
% OUTPUTS:
%  particles    - particles at time t
%  weights      - particle weights at time t
%  mu           - state mean at time t
%  sigma        - state covar at time t 
%  predParticles- predicted particles from PF prediction step
%  predWeights  - predicted particle weights from PF prediction step
%  predMu       - state mean from PF prediction step
%  predSigma    - state covar from PF prediction step
%  pOfZ         - observation likelihood
% You should almost certainly leave the header the way it is.  In particular, you need to compute all return values.

% TODO: Remove this line

% ----------------------------------------------------------------
% ----------------------------------------------------------------
% Prediction step
% ----------------------------------------------------------------
% ----------------------------------------------------------------
predWeights = zeros(numParticles, 1);
n_factor = 0; % normalization factor
for i=1:numParticles
    predParticles(i, :)  = sampleOdometry(u, particles(i, :), alphas);
    zHat_i = observation(predParticles(i, :), FIELDINFO, markerId);
    weight = likelihood(beta, minimizedAngle(z(1) - zHat_i(1)));
    predWeights(i) = weight;
    n_factor = n_factor + weight;
    
end
[predMu, predSigma] = meanAndVariance(transpose(predParticles), numParticles);
pOfZ = mean(predWeights);
predWeights = predWeights / n_factor; % normalize weights

% ----------------------------------------------------------------
% ----------------------------------------------------------------
% Correction step
% ----------------------------------------------------------------
% ----------------------------------------------------------------

[particles, weights] = resample(predParticles, predWeights);

[mu, sigma] = meanAndVariance(transpose(particles), numParticles); 
