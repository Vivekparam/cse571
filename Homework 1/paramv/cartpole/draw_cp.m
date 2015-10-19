%% draw_cp.m
% *Summary:* Draw the cart-pole system with reward, applied force, and 
% predictive uncertainty of the tip of the pendulum
%
%    function draw_cp(state, force, text1, stateSamples)
%
%
% *Input arguments:*
%
%		state           [x,v,dtheta,theta] (D x 1)
%   force           force applied to cart
%   text1           text field 1
%   stateSamples    sample states for display (D x N)
%
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Modified by Arunkumar Byravan for CSE 571: Probabilistic Robotics 
% (Autumn 2015) at  University of Washington

function draw_cp(state, force, text1, stateSamples)
%% Code

% Cartpole params
l = 0.5;
xmin = -5; 
xmax = 5;    
height = 0.1;
width  = 0.3;
maxU = 10;
yoffset = 2.0;

%% Iterate over mean and samples
if (nargin <= 3)
  stateSamples = [];
end
allStates = [state stateSamples];

% Plot a cartpole for each of the states
hold on; % Start holding on the figure
plot([xmin, xmax], [-height-0.03, -height-0.03] + yoffset,'k','linewidth',2); % Plot base line
for k = 1:size(allStates,2)
  
    % Get x and theta
    x = allStates(1,k);
    theta = allStates(4,k);

    % Compute positions 
    cart = [ x + width,  height 
             x + width, -height
             x - width, -height
             x - width,  height
             x + width,  height];
    pendulum = [x, 0; x+2*l*sin(theta), -cos(theta)*2*l];

    % Add offset in y to shift it up
    cart(:,2) = cart(:,2) + yoffset;
    pendulum(:,2) = pendulum(:,2) + yoffset;
    
    % Set transparency for plot objects
    if k == 1; Alpha = 1; else Alpha = 0.25; end;

    % Plot the cart-pole
    p1 = fill(cart(:,1), cart(:,2),'k','edgecolor','k'); alpha(p1, Alpha);
    patchline(pendulum(:,1), pendulum(:,2),'facecolor', 'r',...
      'edgecolor','r','linewidth',4,'edgeAlpha',Alpha,'faceAlpha',Alpha); 

    % Plot the joint and the tip
    dx = [-0.06, 0.06]; dy = [0, 0];
    patchline(dx + x, dy + yoffset,'facecolor', 'y',...
      'edgecolor','y','linewidth',5,'edgeAlpha',Alpha,'faceAlpha',Alpha);
    patchline(dx + pendulum(2,1), dy + pendulum(2,2),'facecolor', 'y',...
      'edgecolor','y','linewidth',5,'edgeAlpha',Alpha,'faceAlpha',Alpha);  
end

%% Common items

% Plot force
plot([0 force/maxU*(xmax-2)],[0.5, 0.5],'g','linewidth',10)
text(0,0.55,'applied force','fontweight','bold','FontName','Times New Roman','fontsize',14);

% Text
if exist('text1','var')
  text(0,0, text1,'fontweight','bold','FontName','Times New Roman','fontsize',14);
end

%% Axis limits
set(gca,'XLim',[xmin xmax],'YLim',[-0.1 3.5]);
set(gca,'DataAspectRatio',[1 1 1]);
%axis square;
axis off;
drawnow;