clear all; close all; clc;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1.0);
set(groot, 'defaultLineMarkerSize', 2);
set(groot, 'defaultFigureUnits', 'centimeters');
% set(groot, 'defaultFigurePosition', [0, 0, 8.5, 6.0]); %single column
set(groot, 'defaultFigurePosition', [100, 100, 16.0, 10.0]); %double column

% Ellipse100_50 = ParticleData('../CSV_Data/Cylinders/Ellipse100_50/file', 400, 'ellipse', 1, 2, 1);
% fig1 = figure; hold on;
% Ellipse100_50.PlotParticlesVorticityRange(fig1, 1, [-5 5], [4 12], [2 6]);
% axis equal;
% set(gca, 'Visible', 'off');
% set(gca, 'FontSize', 11); % Adjust axes font size
% set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
% exportgraphics(gcf, 'LatexFigures/Ellipse100_50.pdf', 'ContentType', 'image', 'Resolution', 300);

% Triangle100_50 = ParticleData('../CSV_Data/Cylinders/Triangle100_50/file', 877, 'triangle', 1, 2, 1);
% fig1 = figure; hold on;
% Triangle100_50.PlotParticlesVorticityRange(fig1, 1, [-5 5], [4 12], [2 6]);
% axis equal;
% set(gca, 'Visible', 'off');
% set(gca, 'FontSize', 11); % Adjust axes font size
% set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
% exportgraphics(gcf, 'LatexFigures/Triangle100_50.pdf', 'ContentType', 'image', 'Resolution', 300);

% Ellipse1000_50 = ParticleData('../CSV_Data/Cylinders/Ellipse1000_50/file', 355, 'ellipse', 1, 2, 1);
% fig1 = figure; hold on;
% Ellipse1000_50.PlotParticlesVorticityRange(fig1, 1, [-10 10], [4 12], [2 6]);
% axis equal;
% set(gca, 'Visible', 'off');
% set(gca, 'FontSize', 11); % Adjust axes font size
% set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
% exportgraphics(gcf, 'LatexFigures/Ellipse1000_50.pdf', 'ContentType', 'image', 'Resolution', 300);

Triangle1000_50 = ParticleData('../CSV_Data/Cylinders/Triangle1000_50/file', 886, 'ellipse', 1, 2, 1);
fig1 = figure; hold on;
Triangle1000_50.PlotParticlesVorticityRange(fig1, 1, [-10 10], [4 12], [2 6]);
axis equal;
set(gca, 'Visible', 'off');
set(gca, 'FontSize', 11); % Adjust axes font size
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
exportgraphics(gcf, 'LatexFigures/Triangle1000_50.pdf', 'ContentType', 'image', 'Resolution', 300);
