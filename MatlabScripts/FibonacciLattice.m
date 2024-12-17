clear all; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 3);
set(groot, 'defaultLineMarkerSize', 2);
% This is for exportgraphics
set(groot, 'defaultFigureUnits', 'centimeters');
% set(groot, 'defaultFigurePosition', [0, 0, 8.5, 6.0]); %single column
set(groot, 'defaultFigurePosition', [100, 100, 16.0, 16.0]); %double column

phi = (1 + sqrt(5)) / 2;
R = 1;
Npoints = 500;
xyz = zeros(Npoints, 3);

figure; hold on;

for i = 1:Npoints
    % theta
    latitude = asin((2 * i - Npoints - 1) / Npoints);
    % phi
    longitude = 2 * pi * i / phi;

    %cartesian coordinates
    x = R * cos(longitude) * cos(latitude);
    y = R * sin(longitude) * cos(latitude);
    z = R * sin(latitude);
    xyz(i, :) = [x, y, z];

    % plot3(xyz(i, 1), xyz(i, 2), xyz(i, 3), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

end

plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 'k.', 'MarkerSize', 20, 'MarkerFaceColor', 'k');

[x y z] = sphere(1000);
x = 1 * x;
y = 1 * y;
z = 1 * z;
surf(x, y, z, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
axis equal; axis tight;
view(3);
set(gca, 'Visible', 'off');
print('LatexFigures/FibonacciLattice', '-dpng', '-r1000');
