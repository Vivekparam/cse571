function [mu, sigma] = meanAndVariance( state, numSamples)
% Assumes equally weighted particles.
mu = mean(state')';

% orientation is a bit more tricky.
sinSum = 0;
cosSum = 0;
for s=1:numSamples,
  cosSum = cosSum + cos( state( 3,s));
  sinSum = sinSum + sin( state( 3,s));
end
mu( 3) = atan2( sinSum, cosSum);

% Compute covariance.
zeroMean = state - repmat(mu,1,numSamples);
for s=1:numSamples,
  zeroMean( 3,s) = minimizedAngle( zeroMean( 3,s));
end

sigma = zeroMean(:,:)*zeroMean(:,:)' / numSamples;
