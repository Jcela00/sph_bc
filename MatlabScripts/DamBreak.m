clear all; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1);
set(groot, 'defaultLineMarkerSize', 2);
set(groot, 'defaultFigureUnits', 'centimeters');
set(groot, 'defaultFigurePosition', [100, 100, 16.0, 16.0]); %double column

rho0 = 1; dim = 2; Hconst = 1;

filenums = [160 180 200 210 220];

for n = 1:length(filenums)
    interact = ParticleData('../CSV_Data/DamBreakInteract/file', filenums(n), ['lolo'], rho0, dim, Hconst);
    nointeract = ParticleData('../CSV_Data/DamBreakNoInteract/file', filenums(n), ['lolo'], rho0, dim, Hconst);

    vely_max = 0.02;
    vely_min = -0.02;

    markersize = 1;
    fig1 = figure; hold on;
    PlotV(interact, fig1, markersize, [vely_min, vely_max]);
    axis equal; axis tight;
    exportgraphics(fig1, ['LatexFigures/DamBreakInteract' num2str(filenums(n)) '.pdf'], 'ContentType', 'vector', 'Resolution', 300);

    fig2 = figure; hold on;
    PlotV(nointeract, fig2, markersize, [vely_min, vely_max]);
    axis equal; axis tight;
    exportgraphics(fig2, ['LatexFigures/DamBreakNoInteract' num2str(filenums(n)) '.pdf'], 'ContentType', 'vector', 'Resolution', 300);

end

function PlotV(data, fig, mrksz, vrange)
    figure(fig);
    custom_colormap = paraview_red_blue_colormap(255);
    vmax = vrange(2);
    vmin = vrange(1);

    for k = 1:data.Npart

        if (data.Type(k) == data.TypeFluid)
            vel_y = data.Velocity(k, 2);
            vel_y = max(vmin, min(vmax, vel_y));

            color_index = round(1 + ((255 - 1) / (vmax - vmin)) * (vel_y - vmin));
            plot(data.Position(k, 1), data.Position(k, 2), 'o', 'MarkerEdgeColor', custom_colormap(color_index, :), 'MarkerFaceColor', custom_colormap(color_index, :), 'MarkerSize', mrksz);
        else
            plot(data.Position(k, 1), data.Position(k, 2), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', mrksz);

        end

    end

    fig.Colormap = custom_colormap;

    cb = colorbar;
    set(cb, 'TickLabelInterpreter', 'latex');
    clim([vrange(1), vrange(2)]);
    % title(obj.PlotName);

end

function cmap = paraview_red_blue_colormap(n)

    white_range = [-0.4, 0.4]; % Default white range

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
