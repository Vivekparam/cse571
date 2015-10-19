function l = likelihood( sigma, innovation)
% Gaussian likelihood
propFactor = 1 / sqrt(det(2 * pi * sigma));
l = propFactor * exp( -0.5 * innovation' * inv(sigma) * innovation);



    
