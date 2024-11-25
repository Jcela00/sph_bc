clear all; close all; clc;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1);
set(groot, 'defaultLineMarkerSize', 8);
% This is for exportgraphics
set(groot, 'defaultFigureUnits', 'centimeters');
% set(groot, 'defaultFigurePosition', [0, 0, 8.5, 6.0]); %single column
set(groot, 'defaultFigurePosition', [0, 0, 16.0, 10.0]); %double column

virtual_color = [0.8, 0.8, 0.8];
virtual_size = 6;
marker_size = 4;
vecsize = 0.3;

Ncircle = 22;
Nflat = 22;
Nsine = 32;
% CIRCLE CURVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Npoints = Ncircle;
R = 1.0;
perimeter = pi * R;
rf = 1.0; % refinement factor
darc = pi * R / Npoints;
dx = rf * darc

points = zeros(Npoints + 1, 2);
normals = zeros(Npoints + 1, 2);
curvature = 1 / R;
Vref = dx ^ 2;
rref = Volume2Radius(Vref);

theta = 0;
dtheta = pi / (Npoints);

for k = 1:Npoints + 1

    points(k, :) = R * [cos(theta), sin(theta)];
    normals(k, :) = [cos(theta), sin(theta)];
    theta = theta + dtheta;
end

fig1 = figure; hold on;
% PlotCircle(fig1, [0, 0], R, 1, '-');
thfine = linspace(0, theta - dtheta, 1000);
plot(R * cos(thfine), R * sin(thfine), '-k');
plot(points(:, 1), points(:, 2), 'ko', 'MarkerSize', marker_size, 'MarkerFaceColor', 'k');
quiver(points(:, 1), points(:, 2), normals(:, 1), normals(:, 2), vecsize, 'k');

for k = 1:Npoints + 1
    [vol] = curvedVol(curvature, dx);
    r = Volume2Radius(vol);
    ratio = r ./ rref;
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

set(gca, 'Visible', 'off');
% set(gcf, 'Units', 'inches');
% screenposition = get(gcf, 'Position');
% set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
% print(gcf, 'LatexFigures/circle.pdf', '-dpdf', '-bestfit');

set(gca, 'FontSize', 11); % Adjust axes font size
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
exportgraphics(gcf, 'LatexFigures/circle.pdf', 'ContentType', 'vector', 'Resolution', 300);

% FLAT WALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Npoints = Nflat;
L = pi * R;
darc = L / (Npoints - 1);

points = zeros(Npoints, 2);
normals = zeros(Npoints, 2);
curvature = 0;
dx = darc
Vref = dx ^ 2;
rref = Volume2Radius(Vref);

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
    ratio = r ./ rref;
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
set(gca, 'Visible', 'off');
set(gca, 'FontSize', 11); % Adjust axes font size
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
% print(gcf, 'LatexFigures/convex.pdf', '-dpdf');
exportgraphics(gcf, 'LatexFigures/flat.pdf', 'ContentType', 'vector', 'Resolution', 300);

% SINE WAVE CURVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Npoints = Nsine;
period = 2 * pi;
L = sqrt(2) * ellipticE(period, 1/2);
darc = L / (Npoints - 1);
dx = darc
Vref = dx ^ 2;
rref = Volume2Radius(Vref);
t = linspace(0, period, Npoints);

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

tfine = linspace(0, period, 1000);
figure; hold on;
plot(tfine, sin(tfine), '-k');
plot(positions(:, 1), positions(:, 2), 'ko', 'MarkerSize', marker_size, 'MarkerFaceColor', 'k');
quiver(positions(:, 1), positions(:, 2), normals(:, 1), normals(:, 2), vecsize, 'k');

for k = 1:Npoints
    vol = curvedVol(curvature(k), darc);
    curvature(k);
    r = Volume2Radius(vol);
    ratio = r ./ rref;
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
set(gca, 'Visible', 'off');
set(gca, 'FontSize', 11); % Adjust axes font size
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
exportgraphics(gcf, 'LatexFigures/sine.pdf', 'ContentType', 'vector', 'Resolution', 300);

% figure;
% plot(positions(:, 1), curvature);

%% check that arcs are equispaced
% arc1 = sqrt(2) * ellipticE(positions(1:end - 1, 1), 1/2);
% arc2 = sqrt(2) * ellipticE(positions(2:end, 1), 1/2);
% darcs = arc2 - arc1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FLAT WITH LABELS
% FLAT WALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
virtual_size = 15;
pbound = [0 0];
normals = [0 1];
pfluid = [1 1.2];
pfluid_proj = [pfluid(1) 0];
dx = 1;
rref = Volume2Radius(dx ^ 2);
fig4 = figure; hold on;
plot(pbound(:, 1), pbound(:, 2), 'ko', 'MarkerSize', virtual_size, 'MarkerFaceColor', 'k');
lfline = [pfluid; pfluid_proj];
rab = [pbound; pfluid];
plot(lfline(:, 1), lfline(:, 2), 'k:');
plot(rab(:, 1), rab(:, 2), '--k');
yline(0, 'k', 'LineWidth', 3);
PlotCircle(fig4, pfluid, 3 * dx, 1, '--');
midpoint_rab = 0.5 * (pbound + pfluid) - [0.1 0];
text(midpoint_rab(1), midpoint_rab(2), '    $\mathbf{r_{ab}}$', 'FontSize', 16, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
midpoint_fluid_proj = 0.5 * (pfluid + pfluid_proj) + [0.1 0];
text(midpoint_fluid_proj(1), midpoint_fluid_proj(2), '    $l_f$', 'FontSize', 16, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
% quiver(pbound(:, 1), pbound(:, 2), normals(:, 1), normals(:, 2), vecsize, 'k');

for k = 1:1
    [vol] = curvedVol(0, 1);
    r = Volume2Radius(vol);
    ratio = r ./ rref;
    first = pbound(k, :) - 0.5 * dx * normals(k, :);
    second = pbound(k, :) - 1.5 * dx * normals(k, :);
    third = pbound(k, :) - 2.5 * dx * normals(k, :);
    sizeofpoint = virtual_size * ratio;

    linemarker2third = [pbound; third];
    plot(linemarker2third(:, 1), linemarker2third(:, 2), 'k');
    midpointMarkerfirst = 0.5 * (pbound + first) + [0.1 0];
    text(midpointMarkerfirst(1), midpointMarkerfirst(2), '    $l_{w1}=\frac{1}{2}\Delta x$', 'FontSize', 16, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    middleMarkersecond = 0.5 * (first + second) + [0.1 0];
    text(middleMarkersecond(1), middleMarkersecond(2), '    $l_{w2}=\frac{3}{2}\Delta x$', 'FontSize', 16, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    middleMarkerthird = 0.5 * (second + third) + [0.1 0];
    text(middleMarkerthird(1), middleMarkerthird(2), '    $l_{w3}=\frac{5}{2}\Delta x$', 'FontSize', 16, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

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

theta = -pi / 4;
circlepoint = pfluid + 3 * dx * [cos(0), sin(0)];
line3h = [pfluid; circlepoint];
midpoint = 0.5 * (pfluid + circlepoint) + [0.0 0.15];
text(midpoint(1), midpoint(2), '    $3h$', 'FontSize', 16, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
plot(line3h(:, 1), line3h(:, 2), '--k');
plot(pfluid(1), pfluid(2), 'ok', 'MarkerSize', virtual_size, 'MarkerFaceColor', 'w');
% axis([-2 4 -4 * dx 4 * dx]);
axis equal;
set(gca, 'Visible', 'off');
set(gca, 'FontSize', 11); % Adjust axes font size
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
exportgraphics(gcf, 'LatexFigures/scheme.pdf', 'ContentType', 'vector', 'Resolution', 300);

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
