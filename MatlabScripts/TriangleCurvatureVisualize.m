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

rho0 = 1; dim = 2; Hconst = 1;

lowK = ParticleData('../CSV_Data/TriangleEquilateral_new_differential_lowKmax_60_40_1rf_1prc/file', 0, ['lolo'], rho0, dim, Hconst);
normalK = ParticleData('../CSV_Data/TriangleEquilateral_new_differential_60_40_1rf_1prc/file', 0, ['lolo'], rho0, dim, Hconst);
highK = ParticleData('../CSV_Data/TriangleEquilateral_new_differential_highKmax_60_40_1rf_1prc/file', 0, ['lolo'], rho0, dim, Hconst);

datasets = {normalK, lowK, highK};
titles = {'Normal', 'Low $K$', 'High $K$'};

for k = 1:length(datasets)
    data = datasets{k};

    obstaclePos = data.Position(data.ObstacleIndexes, 1:2);
    normals = data.Normal(data.ObstacleIndexes, 1:2);
    volumes = data.Volumes(data.ObstacleIndexes, 1:3);

    fluidVol = data.Volumes(data.FluidIndexes, 1);
    dx = fluidVol(1, 1);

    Nparticles = length(obstaclePos(:, 1));

    Vref = dx * dx;
    mrksz = 8;

    fig = figure; hold on;

    for n = 1:Nparticles
        plot(obstaclePos(n, 1), obstaclePos(n, 2), 'ok', 'MarkerSize', mrksz / 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        first = obstaclePos(n, :) - 0.5 * dx * normals(n, :);
        second = obstaclePos(n, :) - 1.5 * dx * normals(n, :);
        third = obstaclePos(n, :) - 2.5 * dx * normals(n, :);

        ratio = volumes(n, :) ./ Vref;

        size = mrksz * ratio;

        % plot line joining marker and third

        if (size(1) > 0)
            plot(first(1), first(2), 'o', 'MarkerSize', size(1), 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
            marker2third = [obstaclePos(n, :); first];
        end

        if (size(2) > 0)
            plot(second(1), second(2), 'o', 'MarkerSize', size(2), 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
            marker2third = [obstaclePos(n, :); second];
        end

        if (size(3) > 0)
            plot(third(1), third(2), 'o', 'MarkerSize', size(3), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
            marker2third = [obstaclePos(n, :); third];
        end

        plot(marker2third(:, 1), marker2third(:, 2), '-k', 'LineWidth', 0.5);

    end

    % title(titles{k});
    axis equal; axis tight;

    set(gca, 'Visible', 'off');
    set(gca, 'FontSize', 11);
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
    exportgraphics(fig, ['LatexFigures/TriangleSketch' num2str(k) '.pdf'], 'ContentType', 'vector', 'Resolution', 300);

end

% filenums = [0 20 40 60 80 100 120 140 160 180 200 220];

% for n = 1:length(filenums)
%     flower = ParticleData('../CSV_Data/Flower/file', filenums(n), ['lolo'], rho0, dim, Hconst);
%     fig = figure; hold on;
%     mrksz = 6;
%     fluidPos = flower.Position(flower.FluidIndexes, 1:2);
%     fluidTag = flower.Volumes(flower.FluidIndexes, 2);

%     obstaclePos = flower.Position(flower.ObstacleIndexes, 1:2);

%     fluid1 = fluidPos(fluidTag == 0, :);
%     fluid2 = fluidPos(fluidTag == 100, :);
%     fluid3 = fluidPos(fluidTag == 200, :);

%     darkblue = [0 0.4470 0.7410];
%     darkerblue = [0 0 153] / 255;
%     orange = [0.8500 0.3250 0.0980];
%     darkyellow = [0.9290 0.6940 0.1250];
%     darkpurple = [0.4940 0.1840 0.5560];
%     mediumgreen = [0.4660 0.6740 0.1880];
%     lightblue = [0.3010 0.7450 0.9330];
%     darkred = [0.6350 0.0780 0.1840];
%     darkerred = [153 0 0] / 255;
%     darkergreen = [0 153 0] / 255;

%     cmap = turbo(256);
%     color1 = cmap(20, :);
%     color2 = cmap(180, :);
%     color3 = cmap(end - 20, :);

%     % Reds:
%     Crimson = [220, 20, 60] / 255;
%     Scarlet = [255, 36, 0] / 255;
%     CherryRed = [255, 28, 28] / 255;
%     ElectricBlue = [125, 249, 255] / 255;
%     RoyalBlue = [65, 105, 225] / 255;
%     Azure = [0, 127, 255] / 255;
%     LimeGreen = [50, 205, 50] / 255;
%     Emerald = [80, 200, 120] / 255;
%     NeonGreen = [57, 255, 20] / 255;

%     color1 = CherryRed;
%     color2 = Azure;
%     color3 = NeonGreen;

%     % plot(obstaclePos(:, 1), obstaclePos(:, 2), 'ok', 'MarkerSize', mrksz, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
%     % plot(fluid1(:, 1), fluid1(:, 2), 'o', 'MarkerSize', mrksz, 'MarkerFaceColor', darkred, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
%     % plot(fluid2(:, 1), fluid2(:, 2), 'o', 'MarkerSize', mrksz, 'MarkerFaceColor', darkyellow, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
%     % plot(fluid3(:, 1), fluid3(:, 2), 'o', 'MarkerSize', mrksz, 'MarkerFaceColor', darkblue, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
%     plot(obstaclePos(:, 1), obstaclePos(:, 2), '.k', 'MarkerSize', mrksz);
%     plot(fluid1(:, 1), fluid1(:, 2), '.', 'MarkerSize', mrksz, 'Color', color1);
%     plot(fluid2(:, 1), fluid2(:, 2), '.', 'MarkerSize', mrksz, 'Color', color2);
%     plot(fluid3(:, 1), fluid3(:, 2), '.', 'MarkerSize', mrksz, 'Color', color3);

%     axis equal; axis tight;
%     set(gca, 'Visible', 'off');
%     set(gca, 'FontSize', 11);
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
%     exportgraphics(fig, ['LatexFigures/Flower' num2str(filenums(n)) '.pdf'], 'ContentType', 'vector', 'Resolution', 300);

% end
