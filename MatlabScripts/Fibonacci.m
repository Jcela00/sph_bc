clear all; clc; close all;
Npoints = 3000;
phi = (1 + sqrt(5)) / 2;
R = 1;
xyz = zeros(Npoints, 3);

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
end

figure;
plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 'or');
axis equal;
