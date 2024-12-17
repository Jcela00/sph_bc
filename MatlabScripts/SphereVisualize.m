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
set(groot, 'defaultFigurePosition', [100, 100, 16.0, 10.0]); %double column

rho0 = 1; dim = 3; Hconst = 1;

sphere = ParticleData('../CSV_Data/Sphere_new_summation_80_40_40_1rf_1prc/file', 3000, ['lolo'], rho0, dim, Hconst);

% Extract particle data
pos = sphere.Position(sphere.FluidIndexes, :);
x = pos(:, 1);
y = pos(:, 2);
z = pos(:, 3);

vel = sphere.Velocity(sphere.FluidIndexes, :);
u = vel(:, 1);
v = vel(:, 2);
w = vel(:, 3);

% dimensions of the box 20x10x10

Nx = 160;
Ny = 80;

% Plot v in the z=5 plane
xq = linspace(0, 20, Nx);
yq = linspace(0, 10, Ny);
zq = 5;
[xq, yq, zq] = meshgrid(xq, yq, zq);
xq = squeeze(xq);
yq = squeeze(yq);
zq = squeeze(zq);

vinterp = griddata(x, y, z, v, xq, yq, zq, 'natural');

levels = [-0.015, -0.01, -0.005, 0.005, 0.01, 0.015];

figure; hold on;
pcolor(xq, yq, vinterp, 'EdgeColor', 'none');
shading interp;
contour(xq, yq, vinterp, levels, 'k');
PlotSphere([10, 5]);
cmap = BlueYellowRedCmap(256);
colormap(cmap);
% c = colorbar;
% set(c, 'TickLabelInterpreter', 'latex');
axis equal; axis tight;
xlabel('$x$'); ylabel('$y$');
title('$v(x, y, z=5)$');
caxis([-0.015, 0.015]);
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/Sphere_v_z.pdf'], 'ContentType', 'vector', 'Resolution', 300);

% Plot w in the y=5 plane
xq = linspace(0, 20, Nx);
yq = 5;
zq = linspace(0, 10, Ny);
[xq, yq, zq] = meshgrid(xq, yq, zq);
xq = squeeze(xq);
yq = squeeze(yq);
zq = squeeze(zq);

winterp = griddata(x, y, z, w, xq, yq, zq, 'natural');
uinterp = griddata(x, y, z, u, xq, yq, zq, 'natural');

figure; hold on;
pcolor(xq, zq, winterp, 'EdgeColor', 'none');
shading interp;
contour(xq, zq, winterp, levels, 'k');
PlotSphere([10, 5]);
cmap = BlueYellowRedCmap(256);
colormap(cmap);
% c = colorbar;
% set(c, 'TickLabelInterpreter', 'latex');
axis equal; axis tight;
xlabel('$x$'); ylabel('$z$');
title('$w(x, y=5, z)$');
caxis([-0.015, 0.015]);
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/Sphere_w_y.pdf'], 'ContentType', 'vector', 'Resolution', 300);

levels_u = [0, 0.03, 0.06, 0.09, 0.10, 0.11, 0.12];
figure; hold on;
pcolor(xq, zq, uinterp, 'EdgeColor', 'none');
shading interp;
contour(xq, zq, uinterp, levels_u, 'k');
PlotSphere([10, 5]);
cmap = BlueGreenYellowRed(256);
colormap(cmap);
% c = colorbar;
% set(c, 'TickLabelInterpreter', 'latex');
axis equal; axis tight;
xlabel('$x$'); ylabel('$z$');
title('$u(x, y=5, z)$');
caxis([0 0.12]);
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/Sphere_u_y.pdf'], 'ContentType', 'vector', 'Resolution', 300);

% Plot u in the x=10 plane
xq = 10;
yq = linspace(0, 10, Ny);
zq = linspace(0, 10, Ny);
[xq, yq, zq] = meshgrid(xq, yq, zq);
xq = squeeze(xq);
yq = squeeze(yq);
zq = squeeze(zq);

uinterp = griddata(x, y, z, u, xq, yq, zq, 'natural');

figure; hold on;
pcolor(yq, zq, uinterp, 'EdgeColor', 'none');
shading interp;
contour(yq, zq, uinterp, levels_u, 'k');
PlotSphere([5, 5]);
cmap = BlueGreenYellowRed(256);
colormap(cmap); c = colorbar;
set(c, 'TickLabelInterpreter', 'latex');
axis equal; axis tight;
xlabel('$y$'); ylabel('$z$');
title('$u(x=10, y, z)$');
caxis([0 0.12]);
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/Sphere_u_x.pdf'], 'ContentType', 'vector', 'Resolution', 300);

% plot v,w at the x=9.5 plane
xq = 9.5;
yq = linspace(0, 10, Ny);
zq = linspace(0, 10, Ny);
[xq, yq, zq] = meshgrid(xq, yq, zq);
xq = squeeze(xq);
yq = squeeze(yq);
zq = squeeze(zq);

uinterp = griddata(x, y, z, u, xq, yq, zq, 'natural');
vinterp = griddata(x, y, z, v, xq, yq, zq, 'natural');
winterp = griddata(x, y, z, w, xq, yq, zq, 'natural');

levels = [-0.015, -0.01, -0.005, 0, 0.005, 0.01, 0.015];

figure; hold on;
pcolor(yq, zq, vinterp, 'EdgeColor', 'none');
shading interp;
contour(yq, zq, vinterp, levels, 'k');
cmap = BlueYellowRedCmap(256);
colormap(cmap); c = colorbar;
set(c, 'TickLabelInterpreter', 'latex');
axis equal; axis tight;
xlabel('$y$'); ylabel('$z$');
title('$v(x=9.5, y, z)$');
caxis([-0.015, 0.015]);
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/Sphere_v_x.pdf'], 'ContentType', 'vector', 'Resolution', 300);

figure; hold on;
pcolor(yq, zq, winterp, 'EdgeColor', 'none');
shading interp;
contour(yq, zq, winterp, levels, 'k');
% PlotSphere([5, 5]);
cmap = BlueYellowRedCmap(256);
colormap(cmap); c = colorbar;
set(c, 'TickLabelInterpreter', 'latex');
axis equal; axis tight;
xlabel('$y$'); ylabel('$z$');
title('$w(x=9.5, y, z)$');
caxis([-0.015, 0.015]);
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/Sphere_w_x.pdf'], 'ContentType', 'vector', 'Resolution', 300);

function [] = PlotSphere(centre)
    sphere_color = [0.7 0.7 0.7];

    cx = centre(1);
    cy = centre(2);
    th = linspace(0, 2 * pi, 100);
    xx = cx + 0.5 * cos(th);
    yy = cy + 0.5 * sin(th);
    fill(xx, yy, sphere_color, 'EdgeColor', 'k', 'LineWidth', 1);
end

function cmap = redwhiteblue(n)
    white_range = [-0.1, 0.1];

    % Define the RGB values for the minimum, central, and maximum colors
    rgb_min = [59, 76, 192] / 255; % Minimum blue
    rgb_mid = [222, 220, 218] / 255; % Central ghost white
    rgb_max = [180, 4, 38] / 255; % Maximum red

    % Define normalized positions for the color interpolation
    x = linspace(-1, 1, n); % Normalized colormap range [-1, 1]
    white_min = white_range(1); % Start of white range
    white_max = white_range(2); % End of white range

    % Initialize the colormap
    cmap = zeros(n, 3);

    % Create the colormap by interpolating RGB values
    for i = 1:3
        % Blue to white
        idx_blue = x <= white_min;
        cmap(idx_blue, i) = interp1([-1, white_min], [rgb_min(i), rgb_mid(i)], x(idx_blue));

        % White region
        idx_white = (x > white_min) & (x < white_max);
        cmap(idx_white, i) = rgb_mid(i);

        % White to red
        idx_red = x >= white_max;
        cmap(idx_red, i) = interp1([white_max, 1], [rgb_mid(i), rgb_max(i)], x(idx_red));
    end

end

function cmap = BlueGreenYellowRed(nColors)
    white = [1, 1, 1]; % White
    light_blue = [0.6, 0.8, 1]; % Light Blue
    green = [0, 0.8, 0]; % Green
    yellow = [1, 1, 0]; % Yellow
    red = [1, 0, 0]; % Red

    cmap = [linspace(white(1), light_blue(1), ceil(nColors / 4))', linspace(white(2), light_blue(2), ceil(nColors / 4))', linspace(white(3), light_blue(3), ceil(nColors / 4))';
            linspace(light_blue(1), green(1), ceil(nColors / 4))', linspace(light_blue(2), green(2), ceil(nColors / 4))', linspace(light_blue(3), green(3), ceil(nColors / 4))';
            linspace(green(1), yellow(1), ceil(nColors / 4))', linspace(green(2), yellow(2), ceil(nColors / 4))', linspace(green(3), yellow(3), ceil(nColors / 4))';
            linspace(yellow(1), red(1), ceil(nColors / 4))', linspace(yellow(2), red(2), ceil(nColors / 4))', linspace(yellow(3), red(3), ceil(nColors / 4))'];
end

function cmap = BlueYellowRedCmap(nColors)
    % Define adjusted colors
    blue = [0, 0, 1]; % Darker blue
    light_blue = [0, 0.8, 1]; % Shiny cyan blue
    white = [1, 1, 1]; % White
    yellow = [1, 0.8, 0]; % Shiny orange
    red = [1, 0, 0];

    % Split evenly between negative and positive sides
    nPart = ceil(nColors / 4);

    % Interpolate colors
    % nPart = floor(nColors / 4);
    blue_part = [linspace(blue(1), light_blue(1), nPart)', linspace(blue(2), light_blue(2), nPart)', linspace(blue(3), light_blue(3), nPart)';
                 linspace(light_blue(1), white(1), nPart)', linspace(light_blue(2), white(2), nPart)', linspace(light_blue(3), white(3), nPart)'];
    red_part = [linspace(white(1), yellow(1), nPart)', linspace(white(2), yellow(2), nPart)', linspace(white(3), yellow(3), nPart)';
                linspace(yellow(1), red(1), nPart)', linspace(yellow(2), red(2), nPart)', linspace(yellow(3), red(3), nPart)'];

    cmap = [blue_part; red_part];

end
