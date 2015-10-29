function [newParticles, newWeights] = resample( particles, weights )
% Given the current particles and weights, sample new particles and weights
% Use the low variance sampler from class.
% INPUTS:
%  particles    - current set of particles
%  weights      - weights for each particle
% 
% OUTPUTS:
% newParticles  - resampled particles (as many as the number of input
%                 particles)
% newWeights    - updated weights for the particles

M = length(particles);
newParticles = zeros(M, 3);

r = rand() * (1 / M);
c = weights(1);
i = 1;

for m = 1:M
    U = r + (m - 1) * (1 / M);
    while U > c
        i = i + 1;
        c = c + weights(i);
    end
    newParticles(m, : ) = particles(i, : );
end

newWeights = ones(M, 1) * (1 / M);
