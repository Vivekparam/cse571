function SLL = gptest_1d(interactive, useMaternKernel, noiseSigma, kernelLengthScale, kernelScaleFactor)
% Input:
%   interactive       - Flag for interactivity (Default: false)
%   useMaternKernel   - Flag for using the matern kernel. By default, this is false
%   noiseSigma        - noise std. deviation (\sigma_n)
%   kernelLengthScale - length scale for the SE kernel and Matern kernel (l)
%   kernelScaleFactor - Scale factor for the SE kernel (\sigma_f)
%                       If the matern kernel is chosen, v = kernelScaleFactor
% TO modify num training points for non-interactive case, set parameter
% "numTrain" in line 46.
if nargin < 1
	interactive = false; % By default, we are in non-interactive mode
end;
if nargin < 2
  useMaternKernel = false; % By default use the Squared Exponential Function
end
if nargin < 3
	noiseSigma = 0.35; %\sigma_n (Noise std. dev)
end;
if nargin < 4
	kernelLengthScale = 1; %\M = kernelLengthScale^-2 (SE), l (Matern)
end;
if nargin < 5
	kernelScaleFactor = 1; %\sigma_f (SE), v (Matern)
end;

% If we use matern kernel, set some default params
if(useMaternKernel)
  kernelLengthScale = 1; % length(l) = 1
  kernelScaleFactor = 1.5; % v = 3/2 (good choices are between 0.5 - 2.5)
end

% same run every time
rng('default');

% the range of x values to show
minX = 0;
maxX = 10;

%% Generate training data if we in non-interactive mode
% Training data is from the sine function with noise added
trainX = []; trainY = [];
if ~interactive
  
  % This is the number of training points
  numTrain = 20;

  % uniform random x values between minX and maxX
  trainX = minX + maxX * rand(1, numTrain);
  
  % The actual data
  theFunction = @(x) sin(x);
  f = theFunction(trainX);
  trainY = zeros(1,numTrain);
  for i=1:numTrain
    trainY(i) = f(i) + noiseSigma * randn(1); % Use SD not variance here
  end

end

%% Generate test queries and ground truth
numTest = 100;
increment=(maxX - minX)/(numTest-1);
testX=minX:increment:maxX; % test queries
if(~interactive)
  testY_gt = theFunction(testX); % groundtruth test data (no-noise added)
end

%% Predict test-results based on the training data
% For the interactive version, you can click on the image and generate
% training points
fg = figure(100); clf;
set(fg, 'Position', [150, 150, 1000, 600]);
if (interactive)
  ylimits = [-5 5];
else
  ylimits = [];
end
while(1)
  
  % If in interactive mode, create a subplot for the ginput
  if(interactive)
    figure(100); 
    subplot(121); hold on; cla;
    xlabel('X'); ylabel('Y'); 
    xlim([minX maxX]); ylim(ylimits); axis square;
    plot(trainX,  trainY, 's', 'MarkerEdgeColor','k', 'MarkerFaceColor','b','MarkerSize',10); % Plot training data	
    title('Click and hit [Enter] add training point. Noise will be added to it.');
    [x,y] = getpts(gca);
    trainX = [trainX x']; %Add to training data
    trainY = [trainY y'+noiseSigma*randn(1,numel(x))]; %Add to training data (noisy signal)
  end
  
  %%% TODO: Condition on the training data and predict test output
  %%% Predict the mean and point-covariance for each test output
  %%% Use the Squared Exponential kernel from the lectures and the
  %%% parameters defined above. 
  % trainX - training inputs, trainY - training outputs
  testY_mean = zeros(1,numTest); % mean prediction for each test point
  testY_var  = ones(1,numTest); % pointwise co-variance
  
  if (useMaternKernel)
    % TODO: Implement conditioning using Matern kernel
    % You might need the functions "besselk" and "gamma"
    % v = kernelScaleFactor, l = kernelLengthScale
    % Use the abs of euclidean distance (r = |x - x'|)
    % Look at eqn 4.14 of: http://www.gaussianprocess.org/gpml/chapters/RW4.pdf

    kernelSize = numel(trainX);
    K = zeros(kernelSize, kernelSize);
    for i=1:kernelSize
        for j=1:kernelSize;
            xi = trainX(i);
            xj = trainX(j);
            K(i, j) = compute_maternk(xi, xj, kernelScaleFactor, kernelLengthScale);
        end
    end
    
    I = eye(kernelSize); % identity matrix
    L = (K + noiseSigma.^2 * I);
    for i=1:numel(testX)
        testXi = testX(i);
        kStar = zeros(1, numel(trainX));
        for j=1:numel(trainX)
           kStar(j) = compute_maternk(trainX(j), testXi, kernelScaleFactor, kernelLengthScale);
        end
        mn = dot(kStar * inv(L), trainY);
        testY_mean(i) = mn;
        variance = compute_maternk(testXi, testXi, kernelScaleFactor, kernelLengthScale) - dot(kStar * inv(L), transpose(kStar));
        testY_var(i) = variance;
    end
  else
    % TODO: Implement the SE kernel predictions
    % kernelLengthScale
    % kernelScaleFactor
    % noiseSigma
    % cov(f(xp), f(xq)) = k(xp, xq) = (-0.5 * magnitude(xp - xq)^2)
    % testX is the test queries. preduct outputs
    kernelSize = numel(trainX);
    K = zeros(kernelSize, kernelSize);
    for i=1:kernelSize
        for j=1:kernelSize;
            xi = trainX(i);
            xj = trainX(j);
            K(i, j) = compute_sek(xi, xj, kernelScaleFactor, kernelLengthScale);
        end
    end
    
    I = eye(kernelSize); % identity matrix
    L = (K + noiseSigma.^2 * I);
    for i=1:numel(testX)
        testXi = testX(i);
        kStar = zeros(1, numel(trainX));
        for j=1:numel(trainX)
           kStar(j) = compute_sek(trainX(j), testXi, kernelScaleFactor, kernelLengthScale);
        end
        mn = dot(kStar * inv(L), trainY);
        testY_mean(i) = mn;
        variance = compute_sek(testXi, testXi, kernelScaleFactor, kernelLengthScale) - dot(kStar * inv(L), transpose(kStar));
        testY_var(i) = variance;
    end
  end
  
  %%% END TODO
  
  %% Display training data
  figure(100);
  subplot(121); hold on; cla;
  xlabel('X'); ylabel('Y'); 
  title('Training data');
  xlim([minX, maxX]);
  if ~isempty(ylimits) ylim(ylimits); end
  axis square;
	plot(trainX,  trainY, 's', 'MarkerEdgeColor','k', 'MarkerFaceColor','b','MarkerSize',10); % Plot training data	
  legend('Noisy training data', 'Location', [0.24, 0.08, 0.1, 0.1]);

  % If we are not in interactive mode, show the ground truth
  if(~interactive)
    plot(testX, testY_gt, 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2); % Plot ground-truth
    legend('Noisy training data', 'Ground truth', 'Location', [0.24, 0.08, 0.1, 0.1]);
  end
  
  %% Display the conditioned result
  figure(100);
  subplot(122); hold on; cla;
  xlabel('X'); ylabel('Y'); 
  title('Test predictions - means and variances');
  xlim([minX maxX]); 
  if ~isempty(ylimits) ylim(ylimits); end
  axis square;

  % Plot results
	plot(trainX,  trainY, 's', 'MarkerEdgeColor','k', 'MarkerFaceColor','b','MarkerSize',10); % Plot training data	
	plot(testX, testY_mean, 'Color', 'k', 'LineWidth', 2); 	% GP mean
  
  % Plot +2*sigma
  plot(testX, testY_mean+2.0*(sqrt(testY_var)), 'Color', 'r', 'LineWidth', 2); 	% GP variance
  plot(testX, testY_mean+2.0*sqrt(testY_var)+noiseSigma, 'Color', 'g', 'LineWidth', 2); % GP variance + noise
  legend('Training data', 'Test mean','Test var','Test var+noise', 'Location', [0.665,0.03,0.15,0.15]);

  % If we are not in interactive mode, show the ground truth
  if(~interactive)
    plot(testX, testY_gt, 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2); % Plot ground-truth
    legend('Training data', 'Test mean','Test var','Test var+noise','Ground truth', 'Location', [0.665,0.03,0.15,0.15]);
  end
  
  % Plot -2*sigma
  plot(testX, testY_mean-2.0*(sqrt(testY_var)), 'Color', 'r', 'LineWidth', 2); 	% GP variance
	plot(testX, testY_mean-2.0*sqrt(testY_var)-noiseSigma, 'Color', 'g', 'LineWidth', 2); % GP variance + noise
  
	%% Compute standardized log loss and return it
  if(~interactive)
    sigmaStarsSquared = testY_var + noiseSigma^2;
    SLL = 0.5 * log(2 * pi * sigmaStarsSquared) + (testY_gt - testY_mean).^2 ./ (2 * sigmaStarsSquared);
    disp(['Mean Standardized Log Loss:', num2str(mean(SLL))]);
    break;
  else
    SLL = [];
  end
  
end

end

