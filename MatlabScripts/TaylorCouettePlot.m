clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 6);
% This is for exportgraphics
set(groot, 'defaultFigureUnits', 'centimeters');
% set(groot, 'defaultFigurePosition', [0, 0, 8.5, 6.0]); %single column
set(groot, 'defaultFigurePosition', [100, 100, 16.0, 10.0]); %double column

rho0 = 1;
dim = 2;
Hconst = 1;

TC = ParticleData(['../CSV_Data/TaylorCouette/TaylorCouette_new_differential_constV_10_10_1rf_1prc/file'], 300, ['20 new'], rho0, dim, Hconst);

mrksz = 4;

figure1 = figure; hold on;
TC.PlotParticles(figure1, mrksz);
xlabel('$x$'); ylabel('$y$');
axis equal; axis tight;
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(figure1, ['LatexFigures/TaylorCouetteParticles.pdf'], 'ContentType', 'vector', 'Resolution', 300);

pos = TC.Position(TC.FluidIndexes, 1:2);
r = sqrt(pos(:, 1) .^ 2 + pos(:, 2) .^ 2);

vel = TC.Velocity(TC.FluidIndexes, 1:2);
params = [0.5 1 1 1];

[ur, uth] = TC_Analytical(r, params);

rfine = linspace(0.5, 1, 1000);
[urfine, uthfine] = TC_Analytical(rfine, params);

figure2 = figure; hold on;
plot(rfine, uthfine, 'k', 'DisplayName', 'Analytical');
plot(r, uth, 'xb', 'DisplayName', 'SPH $h=\frac{1}{20}$');
xlabel('$r$'); ylabel('$u_{\theta}$');
legend('Location', 'best');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(figure2, ['LatexFigures/TaylorCouettePlot.pdf'], 'ContentType', 'vector', 'Resolution', 300);

function [ur, uth] = TC_Analytical(r, params)
    Rin = params(1);
    Rout = params(2);
    Win = params(3);
    Wout = params(4);
    a =- ((Rout ^ 2 * Rin ^ 2) / (Rout ^ 2 - Rin ^ 2)) * (Wout - Win);
    b = (Wout * Rout ^ 2 - Win * Rin ^ 2) / (Rout ^ 2 - Rin ^ 2);
    ur = zeros(length(r), 1);
    uth = a ./ r + b * r;

end
