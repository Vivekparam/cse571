%-------------------------------------------------------
% samples from a multivariate gaussian with given mean and covariance
%-------------------------------------------------------

function s = sample(mean, cov, dimensions)

s=mean+chol(cov)*randn(dimensions,1);
