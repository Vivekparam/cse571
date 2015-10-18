clear all;
close all;
addpath('cartpole/');
settings_cp; % Sets parameters for the simulation, start states etc

%% SETUP
% Cartpole state:
%       x:        [m]     position of cart
%       dx:       [m/s]   velocity of cart
%       dtheta:   [rad/s] angular velocity
%       theta:    [rad]   angle
% Control:
%       u:        [N]     force on cart
% Dimensions of state and control
sdim = 4; udim = 1;

% Dynamics model input:
%       [x,dx,dtheta,sin(theta),cos(theta),u]
% Dynamics prediction: (separate GP for each output)
%       [dx,d^2x,d^2theta,dtheta] % delta-state, not next state
dyni = 6; dyno = 4;

% GP parameters (one per GP) (Try the effect of different parameters)
% 6 dimensions per GP, 4 GP in total
kernelLengthScale = [203.0197  242.9594  240.5070  218.0256
                     33.0219  176.8396  175.9314  178.0185
                     34.6307    7.3903    7.4687   13.0914
                      2.3903    1.0499    0.8433    1.2963
                     31.2894    0.9858    0.7810    1.7216
                    219.1850   24.6355   23.1603   49.9782]; % separate length for each dimension
kernelScaleFactor = [1.1478    1.3658    3.5236    0.7204]; % \sigma_f
noiseSigma        = [0.0143    0.0165    0.0431    0.0145]; % \sigma_n

% Init
numTrainEpochs     = 25; 
numDataPtsPerEpoch = H; % Number of traj points (see settings_cp.m)
numTrajSamples     = 10; 

% Generate an initial rollout on the cartpole using full dynamics
% xx = [x_t, dx_t, dtheta_t, theta_t, sin(theta_t), cos(theta_t), u_t];
% yy = [x_(t+1), dx_(t+1), dtheta_(t+1), theta_(t+1)] (dynamics predictions for training data)
[xx_fulldyn, yy_fulldyn] = rollout(gaussian(mu0, S0), policy, H, plant);
rollout_fulldyn = xx_fulldyn(:, plant.odei)'; % odei = state variable indices
controlTraj     = xx_fulldyn(:, end-udim+1:end)'; % controls

% Use this rollout for initial training data
% NOTE: WE are predicting delta-state, so we use: state_(t+1) - state_t
trainX = xx_fulldyn'; trainX(plant.angi,:) = []; % Remove the angle co-ordinate
trainY = (yy_fulldyn - xx_fulldyn(:,plant.odei))'; 

% Run the training loop. Use the full dynamics (cartpole_dynamics.m) to 
% generate the training dataset at each epoch. 
% Compare known dynamics against predictions from the GP model
for k = 1:numTrainEpochs
  
  % Switch between random and learned policy
  if (~mod(k,4)) % every 4th iteration
    disp('Using learned policy');
    policy = dataset.policy; % Use previously learned policy
    S0 = diag([0.05, 0.01, 0.01, 0.05].^2); % Change initial state distribution
  else
    disp('Using random policy');
    policy = []; policy.maxU = 10; % random policy
    S0 = diag([0.1 0.01 0.01 pi].^2);   % initial state covariance
  end
  
  % Generate a initial state by sampling from a gaussian
  initState = gaussian(mu0, S0); % start state (s_0)

  %% EVALUATE THE GP MODEL AND FULL DYNAMICS USING THE CONTROL TRAJECTORY
  % 1) Get the result from the full dynamics
  % THIS WILL BE USED AS THE TRAINING DATA AT THE END OF THE CURRENT EPOCH
  [xx_fulldyn, yy_fulldyn] = rollout(initState, policy, H, plant);
  rollout_fulldyn = xx_fulldyn(:, plant.odei)'; % odei = state variable indices
  controlTraj     = xx_fulldyn(:, end-udim+1:end)'; % controls / sequence of controls (u0, u1, u2 ... u H-1)  
  
  %% STUDENT TODO:
  % 2a) TODO: Get the result from the GP model
  % Go sequentially along the trajectory
  % Condition on the training data [trainX, trainY] to generate predictions and variances
  % You need to have one GP per output
  % All the 4 variables below need to be populated
  pred_gp_mean    = zeros(dyno, numDataPtsPerEpoch); % mean predictions from the GP [dx, d^2x, d^2theta, dtheta]
  pred_gp_var     = ones(dyno, numDataPtsPerEpoch); % variances for the GP predictions
  rollout_gp      = zeros(dyno, numDataPtsPerEpoch); % state of the system [x, dx, dtheta, theta]
  rollout_gp_var  = ones(dyno, numDataPtsPerEpoch); % variances for the state
  
  % Create the Squared Exponential Kernel
  kernelSize = numel(trainY(1, :));
  K = zeros(4, kernelSize, kernelSize);
  for k=1:4
      for i=1:kernelSize
          for j=1:kernelSize;
              xi = trainX(:, i);
              xj = trainX(:, j);
              K(k, i, j) = compute_sek(xi, xj, kernelScaleFactor(k), kernelLengthScale(:, k));
          end
      end
  end
    
  for i=1:numDataPtsPerEpoch
    testXi =    [rollout_fulldyn(1:3, i)
                sin(rollout_fulldyn(4, i))
                cos(rollout_fulldyn(4, i))
                controlTraj(i)];
    
    % iterate through each GP, getting mean/variance
    % dx, d^2x, d^2theta, dtheta    
    for t=1:dyno
        I = eye(kernelSize); % identity matrix
        L = (squeeze(K(t, :, :)) + noiseSigma(t).^2 * I);
        kStar = zeros(1, numel(trainX(1, :)));
        for j=1:numel(kStar)
           kStar(j) = compute_sek(trainX(:, j), testXi, kernelScaleFactor(t), kernelLengthScale(:, t));
        end
        dmn = dot(kStar * inv(L), trainY(t, :));
        pred_gp_mean(t, i) = dmn;
        rollout_gp(t,i) = rollout_gp(t,i) + dmn;
        dvariance = compute_sek(testXi, testXi, kernelScaleFactor(t), kernelLengthScale(:, t)) - dot(kStar * inv(L), transpose(kStar));
        pred_gp_var(t, i) = dvariance;
        rollout_gp_var(t, i) = rollout_gp_var(t, i) + dvariance;
    end
  end
  
  % 2b) TODO: Sample "numTrajSamples" trajectories using the previously computed means 
  % and variances. These will be used to display how the uncertainty
  % accumulates as we go through a long-term rollout over many steps
  %
  % All trajectories start at initState. 
  % At each prediction step, instead of using the mean-prediction, sample
  % from gaussian around the mean-prediction using the computed variance.
  % Use this sample for the next prediction and repeat.
  % The control trajectory is the same for all trajectories.
  % (initState,control) ->  delta-prediction & variance ->
  % gaussian(meanstate, variance) -> (newState, control) -> .....
  %
  % The sampling from the gaussian is the only extra step compared to 2a)
  % All the 4 variables below need to be populated
  pred_gp_mean_trajs    = zeros(dyno, numDataPtsPerEpoch, numTrajSamples);
  pred_gp_var_trajs     = ones(dyno, numDataPtsPerEpoch, numTrajSamples);
  rollout_gp_trajs      = zeros(dyno, numDataPtsPerEpoch, numTrajSamples); 
  rollout_gp_var_trajs  = ones(dyno, numDataPtsPerEpoch, numTrajSamples);
  
  for i=1:numDataPtsPerEpoch
    testXi =    [rollout_fulldyn(1:3, i)
                sin(rollout_fulldyn(4, i))
                cos(rollout_fulldyn(4, i))
                controlTraj(i)];
    
    % iterate through each GP, getting mean/variance
    % dx, d^2x, d^2theta, dtheta    
    for t=1:dyno
        I = eye(kernelSize); % identity matrix
        L = (squeeze(K(t, :, :)) + noiseSigma(t).^2 * I);
        kStar = zeros(1, numel(trainX(1, :)));
        for j=1:numel(kStar)
           kStar(j) = compute_sek(trainX(:, j), testXi, kernelScaleFactor(t), kernelLengthScale(:, t));
        end
        
        % calc mean and variance
        dmn = dot(kStar * inv(L), trainY(t, :));
        dvariance = compute_sek(testXi, testXi, kernelScaleFactor(t), kernelLengthScale(:, t)) - dot(kStar * inv(L), transpose(kStar));

        % Sample means and variances
        for j=1:numTrajSamples
          sampleMean = gaussian(dmn, dvariance);
          
          pred_gp_mean(t, i) = dmn;
          if(j == 1)
            rollout_gp(t, i, j) = rollout_gp(t,i) + dmn;
          else
              
          end
          % Sample variances
          pred_gp_var(t, i) = dvariance;
          rollout_gp_var(t, i) = rollout_gp_var(t, i) + dvariance;
        end
    end
  end
  
  
  %% END TODO
  % 3) Add the full dynamics predictions to the training data
  % NOTE: WE are predicting delta-state, so we use: state_(t+1) - state_t
  augX = xx_fulldyn'; augX(plant.angi,:) = []; % Remove the angle co-ordinate
  augY = (yy_fulldyn - xx_fulldyn(:,plant.odei))';
  trainX = [trainX augX];
  trainY = [trainY augY];
    
  % 4) Display the predictions from the full dynamics and GP side-by-side
  figure(fg);
  for j= 1:numDataPtsPerEpoch
    
    % Display full dynamics result
    subplot(pltHdls(1)); cla;
    title('Full-Dynamics prediction','fontweight','bold','FontName','Times New Roman','fontsize',14);
    draw_cp(rollout_fulldyn(:,j), controlTraj(:,j), ['T=' num2str(j*dt) 'sec']);
    
    % Display GP result
    subplot(pltHdls(2)); cla;
    title(sprintf('GP-Dynamics prediction \n Training epoch: %d, Num Training data: %d',k,size(trainX,2)-numDataPtsPerEpoch), ...
      'fontweight','bold','FontName','Times New Roman','fontsize',14);
    draw_cp(rollout_gp(:,j), controlTraj(:,j), ['T=' num2str(j*dt) 'sec'], squeeze(rollout_gp_trajs(:,j,:)));
    
    % Show plots of GP prediction
    for i = 1:dyno
      subplot(pltHdls(i+2)); cla; hold on;
      deltaTs = dt:dt:dt*j;
      plot(deltaTs, augY(i,1:j),'g-','Linewidth',2.5);
      plot(deltaTs, pred_gp_mean(i,1:j),'r-','Linewidth',2.5);
      ebar = errorbar(deltaTs, pred_gp_mean(i,1:j), 3*sqrt(pred_gp_var(i,1:j)));
      set(ebar, 'Color', 'r','Linewidth',1.5);
      xlabel('Time (s)','fontweight','bold','fontsize',14,'FontName','Times New Roman');
      ylabel(predTxt{i},'fontweight','bold','fontsize',14,'FontName','Times New Roman');
      xlim([0,totT]);
    end
    
    % Pause to display
    drawnow;
    pause(0.1*dt);
  end
  
  % 5) Compare the mean-predictions from the GP and the full dynamics
  % Leave out the starting point
  rollout_gp(4,:) = mod(rollout_gp(4,:),2*pi);
  rollout_fulldyn(4,:) = mod(rollout_fulldyn(4,:),2*pi);
  sigmaStarsSquared = bsxfun(@plus, rollout_gp_var, (noiseSigma').^2);
  SLL = 0.5 * log(2 * pi * sigmaStarsSquared) + (rollout_fulldyn - rollout_gp).^2 ./ (2 * sigmaStarsSquared);
  disp(['Epoch: ', num2str(k), '. Mean Standardized Log Loss for each GP: ']);
  disp(num2str(mean(SLL,2)'));
  
end