function plotLine(center, angle, length, color)

WAS_HOLD = ishold;

if ~WAS_HOLD
    hold on
end

plot([center(1) center(1)+cos(angle)*length], [center(2) center(2)+sin(angle)*length], 'Color', color);

if ~WAS_HOLD
    hold off
end
