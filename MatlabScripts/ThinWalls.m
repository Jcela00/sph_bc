clear all; close all; clc;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1.0);
set(groot, 'defaultLineMarkerSize', 6);
set(groot, 'defaultFigureUnits', 'centimeters');
% set(groot, 'defaultFigurePosition', [0, 0, 8.5, 6.0]); %single column
set(groot, 'defaultFigurePosition', [100, 100, 16.0, 10.0]); %double column

dx = 1;
vector_size = 0.4;
vector_thickness = 2.0;
random_disp = 0.4;

fig = figure; hold on;

imax = 10;

for i = 1:imax

    for j = 1:4

        y = j * dx;
        x = (i - 1 - 0.5) * dx;

        if ((x ~= 3.5) || (y ~= 1))
            x = x + rand() * random_disp;
            y = y + rand() * random_disp;
            plot(x, y, 'bo', 'MarkerFaceColor', 'b');

            if (isInsideCircle(x, y, 3.5, 1, 3.0 * dx))
                plot(x, y, 'ks', 'MarkerSize', 12);
            end

        else
            plot(x, y, 'bo', 'MarkerFaceColor', 'b');
            plot(x, y, 'ro', 'MarkerSize', 12);

        end

        plot(x, y, 'bo', 'MarkerFaceColor', 'b');

        if (y < 3.5)
            quiver(x, y, 0, 1, vector_size, 'k', 'LineWidth', vector_thickness);
        end

    end

end

yline(0.5 * dx, 'k', 'LineWidth', 2.0);
yline(-0.5 * dx, 'k', 'LineWidth', 2.0);

for i = 1:imax

    for j = 1:4

        y = (-j) * dx;
        x = (i - 1 - 0.5) * dx;

        x = x + rand() * random_disp;
        y = y + rand() * random_disp;

        plot(x, y, 'bo', 'MarkerFaceColor', 'b');

        if (y > -3.5)
            quiver(x, y, 0, -1, vector_size, 'k', 'LineWidth', vector_thickness);
        end

    end

end

for i = 1:imax
    x = dx * (i - 1);

    plot(x, 0.5 * dx, 'ko', 'MarkerFaceColor', 'k');
    quiver(x, 0.5 * dx, 0, 1, vector_size, 'k', 'LineWidth', vector_thickness);

    if (isInsideCircle(x, 0.5 * dx, 3.5, 1, 3.0 * dx))
        plot(x, 0.5 * dx, 'ks', 'MarkerSize', 12);
    end

    plot(x, -0.5 * dx, 'ko', 'MarkerFaceColor', 'k');
    quiver(x, -0.5 * dx, 0, -1, vector_size, 'k', 'LineWidth', vector_thickness);

end

% choosen fluid particle [3.5 1]
radius = 3.0 * dx;

plotCircle(fig, 3.5, 1, radius);

axis equal;
axis tight;
set(gca, 'Visible', 'off');
exportgraphics(gcf, 'LatexFigures/ThinWall.pdf', 'ContentType', 'vector', 'Resolution', 300);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig2 = figure; hold on;

imax = 10;
jmax = 5;

% choosen particle 3.5 1.5
plotCircle(fig2, 3.5, 1.5, 3.0 * dx);

for i = 1:imax

    for j = 1:jmax

        y = (j - 0.5) * dx;
        x = (i - 1 - 0.5) * dx;

        if ((x ~= 3.5) || (y ~= 1.5))
            x = x + rand() * random_disp;
            y = y + rand() * random_disp;
            plot(x, y, 'bo', 'MarkerFaceColor', 'b');

            if (isInsideCircle(x, y, 3.5, 1.5, 3.0 * dx))
                plot(x, y, 'ks', 'MarkerSize', 12);
            end

        else
            plot(x, y, 'bo', 'MarkerFaceColor', 'b');
            plot(x, y, 'ro', 'MarkerSize', 12);

        end

        if (y < 2.5 * dx)
            quiver(x, y, 0, 1, vector_size, 'k', 'LineWidth', vector_thickness);
        else
            quiver(x, y, 0, -1, vector_size, 'k', 'LineWidth', vector_thickness);
        end

    end

end

for i = 1:imax
    x = dx * (i - 1);

    plot(x, 0, 'ko', 'MarkerFaceColor', 'k');
    quiver(x, 0, 0, 1, vector_size, 'k', 'LineWidth', vector_thickness);

    plot(x, jmax * dx, 'ko', 'MarkerFaceColor', 'k');
    quiver(x, jmax * dx, 0, -1, vector_size, 'k', 'LineWidth', vector_thickness);

end

yline(0, 'k', 'LineWidth', 2.0);
yline(jmax * dx, 'k', 'LineWidth', 2.0);

axis equal;
axis tight;
set(gca, 'Visible', 'off');
exportgraphics(gcf, 'LatexFigures/ThinDuct.pdf', 'ContentType', 'vector', 'Resolution', 300);

function plotCircle(fighandle, cx, cy, r)
    % cx, cy: center of the circle
    % r: radius of the circle
    theta = linspace(0, 2 * pi, 100); % Generate 100 points from 0 to 2*pi
    x = cx + r * cos(theta); % X-coordinates of the circle
    y = cy + r * sin(theta); % Y-coordinates of the circle
    figure(fighandle); % Set the current figure to fighandle
    plot(x, y, 'k-', 'LineWidth', 1); % Plot the circle
end

function [isInside] = isInsideCircle(x, y, cx, cy, r)
    % x, y: coordinates of the point
    % cx, cy: center of the circle
    % r: radius of the circle
    isInside = (x - cx) ^ 2 + (y - cy) ^ 2 < r ^ 2;
end
