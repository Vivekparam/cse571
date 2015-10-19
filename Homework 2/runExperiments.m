% Initial settings 
numSteps = 200;
usePf = false; % change this to use particle filter
numParticles = 100;

% Run different experiments
% factors for alphas or betas 
% (note that these factors must themselves be squared with the arrayfun as alphas and betas are in squared variance units)
powers = -3:1:3; % powers of 2 (1/64 to 64)
factorList = arrayfun(@(x) (2^x)^2, powers);
%factorList = factorList(end);
for i = 1:numel(factorList)
	f = factorList(i);
    [mu, sigma, positionError, mahalError, pOfZ] = run(numSteps, usePf, 0.001, true, false, [f,f,f,f], numParticles);
    positionErrors(i) = positionError;
    mahalErrors(i) = mahalError;
    ANEES(i) = mahalError / 3;
    pOfZs(i) = pOfZ;
end

% Flag for particle filter
if usePf
  baseString = 'PF ';
else
  baseString = 'EKF ';
end

% Print results
positionErrors
ANEES
pOfZs

% Generate plots. You might want to change the names of these plots to
% something informative
figure(1); plot(powers, positionErrors); title ([baseString, 'positionErrors']); print('figure1.png', '-dpng');
figure(2); plot(powers, ANEES); title ([baseString, 'ANEES']); print('figure2.png', '-dpng');
figure(3); plot(powers, pOfZs); title ([baseString, 'pOfZs']); print('figure3.png', '-dpng');
%disp('paused...');
%pause();