function [mask] = analyze_background(center, size, disc_radius)
    radius = center(3);
    [col row] = meshgrid(1:size(2), 1:size(1));
    circlePx = (row - center(2)).^2 + ...
        (col - center(1)).^2 <= radius.^2;
    circlePx2 = (row - center(2)).^2 + ...
        (col - center(1)).^2 <= (radius + disc_radius).^2;
    mask = circlePx2 - circlePx;
end