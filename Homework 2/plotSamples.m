function plotSamples(samples, color)

WAS_HOLD = ishold;

if ~WAS_HOLD
    hold on
end

plot( samples(1,:), samples(2,:), '.', 'linewidth', 1, 'Color', color);

if ~WAS_HOLD
    hold off
end
