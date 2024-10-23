clear all; close all; clc;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 600 600]);

virtual_color = [0.8, 0.8, 0.8];
virtual_size = 8;
marker_size = 4;
vecsize = 0.3;

% circle curve
Npoints = 40;
R = 1.0;
perimeter = 2 * pi * R;
rf = 1.0; % refinement factor
darc = 2 * pi * R / Npoints;
dx = rf * darc;

points = zeros(Npoints, 2);
normals = zeros(Npoints, 2);
curvature = 1 / R;

theta = 0;
dtheta = 2 * pi / (Npoints);

for k = 1:Npoints

    points(k, :) = R * [cos(theta), sin(theta)];
    normals(k, :) = [cos(theta), sin(theta)];
    theta = theta + dtheta;
end

fig1 = figure; hold on;
PlotCircle(fig1, [0, 0], R, 1, '-');
plot(points(:, 1), points(:, 2), 'ko', 'MarkerSize', marker_size, 'MarkerFaceColor', 'k');
quiver(points(:, 1), points(:, 2), normals(:, 1), normals(:, 2), vecsize, 'k');

for k = 1:Npoints
    [vol] = curvedVol(curvature, dx);
    r = Volume2Radius(vol);
    ratio = r ./ r(1);
    first = points(k, :) - 0.5 * dx * normals(k, :);
    second = points(k, :) - 1.5 * dx * normals(k, :);
    third = points(k, :) - 2.5 * dx * normals(k, :);
    sizeofpoint = virtual_size * ratio;

    if sizeofpoint(1) > 0
        plot(first(1), first(2), 'ko', 'MarkerSize', sizeofpoint(1), 'MarkerFaceColor', virtual_color);
    end

    if sizeofpoint(2) > 0
        plot(second(1), second(2), 'ko', 'MarkerSize', sizeofpoint(2), 'MarkerFaceColor', virtual_color);
    end

    if sizeofpoint(3) > 0
        plot(third(1), third(2), 'ko', 'MarkerSize', sizeofpoint(3), 'MarkerFaceColor', virtual_color);
    end

end

axis equal;

% FLAT WALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Npoints = 10;
L = Npoints - 1;
darc = L / (Npoints - 1);

points = zeros(Npoints, 2);
normals = zeros(Npoints, 2);
curvature = 0

dx = darc;

for k = 1:Npoints

    points(k, :) = [darc * (k - 1), 0];
    normals(k, :) = [0 1];
end

fig2 = figure; hold on;
plot(points(:, 1), points(:, 2), 'ko', 'MarkerSize', marker_size, 'MarkerFaceColor', 'k');
yline(0, 'k', 'LineWidth', 1);
quiver(points(:, 1), points(:, 2), normals(:, 1), normals(:, 2), vecsize, 'k');

for k = 1:Npoints
    [vol] = curvedVol(curvature, dx);
    r = Volume2Radius(vol);
    ratio = r ./ r(1);
    first = points(k, :) - 0.5 * dx * normals(k, :);
    second = points(k, :) - 1.5 * dx * normals(k, :);
    third = points(k, :) - 2.5 * dx * normals(k, :);
    sizeofpoint = virtual_size * ratio;

    if sizeofpoint(1) > 0
        plot(first(1), first(2), 'ko', 'MarkerSize', sizeofpoint(1), 'MarkerFaceColor', virtual_color);
    end

    if sizeofpoint(2) > 0
        plot(second(1), second(2), 'ko', 'MarkerSize', sizeofpoint(2), 'MarkerFaceColor', virtual_color);
    end

    if sizeofpoint(3) > 0
        plot(third(1), third(2), 'ko', 'MarkerSize', sizeofpoint(3), 'MarkerFaceColor', virtual_color);
    end

end

axis equal;

% sine wave curve
Npoints = 30;
L = sqrt(2) * ellipticE(2 * pi, 1/2);
darc = L / (Npoints - 1);

t = linspace(0, 2 * pi, Npoints);
accumulated_arc = 0;
Nplaced = 0;

positions = zeros(Npoints, 2);
normals = zeros(Npoints, 2);
curvature = zeros(Npoints, 1);
arclength = zeros(Npoints, 1);

% place first particle at t=0
positions(1, :) = [0, 0];
normals(1, :) = [-cos(0), 1];
normals(1, :) = normals(1, :) ./ norm(normals(1, :));
arclength(1) = darc;
curvature(1) = sincurvature(0);

for k = 1:Npoints - 1
    fun = @(phi) paramfun(phi, k, darc);
    phi = fzero(fun, t(k)); % equispaced angles as initial guess

    positions(k + 1, :) = [phi, sin(phi)];
    normals(k + 1, :) = [-cos(phi), 1];
    normals(k + 1, :) = normals(k + 1, :) ./ norm(normals(k + 1, :));
    arclength(k + 1) = darc;
    curvature(k + 1) = sincurvature(phi);
end

tfine = linspace(0, 2 * pi, 1000);
figure; hold on;
plot(tfine, sin(tfine), '-k');
plot(positions(:, 1), positions(:, 2), 'ko', 'MarkerSize', marker_size, 'MarkerFaceColor', 'k');
quiver(positions(:, 1), positions(:, 2), normals(:, 1), normals(:, 2), vecsize, 'k');

for k = 1:Npoints
    vol = curvedVol(curvature(k), darc);
    r = Volume2Radius(vol);
    ratio = r ./ r(1);
    first = positions(k, :) - 0.5 * darc * normals(k, :);
    second = positions(k, :) - 1.5 * darc * normals(k, :);
    third = positions(k, :) - 2.5 * darc * normals(k, :);
    sizeofpoint = virtual_size * ratio;

    if sizeofpoint(1) > 0
        plot(first(1), first(2), 'ko', 'MarkerSize', sizeofpoint(1), 'MarkerFaceColor', virtual_color);
    end

    if sizeofpoint(2) > 0
        plot(second(1), second(2), 'ko', 'MarkerSize', sizeofpoint(2), 'MarkerFaceColor', virtual_color);
    end

    if sizeofpoint(3) > 0
        plot(third(1), third(2), 'ko', 'MarkerSize', sizeofpoint(3), 'MarkerFaceColor', virtual_color);
    end

end

axis equal;

% figure;
% plot(positions(:, 1), curvature);

%% check that arcs are equispaced
% arc1 = sqrt(2) * ellipticE(positions(1:end - 1, 1), 1/2);
% arc2 = sqrt(2) * ellipticE(positions(2:end, 1), 1/2);
% darcs = arc2 - arc1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%

function F = paramfun(phi, i, ds)
    F = sqrt(2) * ellipticE(phi, 1/2) - i * ds;
end

function kappa = sincurvature(phi)
    dy2 = -sin(phi);
    dy = cos(phi);

    kappa = -dy2 / ((1 + dy ^ 2) ^ (3/2));
end

function W = kernel(r, h)

    q = r / h;

    K = 7.0/478.0 / pi / h / h; % 2d normalization

    tmp3 = (3.0 - q) .* (3.0 - q) .* (3.0 - q) .* (3.0 - q) .* (3.0 - q);
    tmp2 = (2.0 - q) .* (2.0 - q) .* (2.0 - q) .* (2.0 - q) .* (2.0 - q);
    tmp1 = (1.0 - q) .* (1.0 - q) .* (1.0 - q) .* (1.0 - q) .* (1.0 - q);

    if (q >= 0.0 & q < 1.0)
        W = K * (tmp3 - 6.0 * tmp2 + 15.0 * tmp1);
    elseif (q >= 1.0 & q < 2.0)
        W = K * (tmp3 - 6.0 * tmp2);
    elseif (q >= 2.0 & q < 3.0)
        W = K * tmp3;
    else
        W = 0.0;
    end

end

function W = kernelGradient(r, h)

    q = r / h;

    K = -5.0 * 7.0/478.0 / pi / h / h / h; % 2d normalization

    tmp3 = (3.0 - q) .* (3.0 - q) .* (3.0 - q) .* (3.0 - q);
    tmp2 = (2.0 - q) .* (2.0 - q) .* (2.0 - q) .* (2.0 - q);
    tmp1 = (1.0 - q) .* (1.0 - q) .* (1.0 - q) .* (1.0 - q);

    if (q >= 0.0 & q < 1.0)
        W = K * (tmp3 - 6.0 * tmp2 + 15.0 * tmp1);
    elseif (q >= 1.0 & q < 2.0)
        W = K * (tmp3 - 6.0 * tmp2);
    elseif (q >= 2.0 & q < 3.0)
        W = K * tmp3;
    else
        W = 0.0;
    end

end

function [vol] = curvedVol(kappa, dx)
    ncol = size(kappa, 2);
    vol = zeros(3, ncol);

    for n = 1:3
        vol(n, :) = 0.5 * (2 * dx + dx * dx * kappa - 2 * n * dx * dx * kappa) * dx;
    end

    vol = max(0, vol);

end

function R = Volume2Radius(Volume)
    R = sqrt(Volume / pi);
end

function PlotCircle(fig, centre, radius, thickness, style)
    Nsamples = 1000;
    theta = linspace(0, 2 * pi, Nsamples);
    x = centre(1) + radius * cos(theta);
    y = centre(2) + radius * sin(theta);
    figure(fig);
    plot(x, y, 'k', 'LineWidth', thickness, 'LineStyle', style);
end

function PlotCircularSector(fig, centre, radius, theta0, thetaN, thickness, style, color)
    Nsamples = 1000;
    theta = linspace(theta0, thetaN, Nsamples);
    x = centre(1) + radius * cos(theta);
    y = centre(2) + radius * sin(theta);
    figure(fig);
    plot(x, y, 'LineWidth', thickness, 'LineStyle', style, 'Color', color);
end
