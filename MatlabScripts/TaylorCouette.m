clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 600 600]);
dp = 0.05;
Nfluid = [80 80 1];
Nboundary = [0, 0, 0];
Hconst = 1;
dim = 2;
rho0 = 1;

TC = ParticleData('../CSV_Data/TaylorCouette_NEW_BC_Summation_80_80_4prc/file', 1000, ['TC 80x80'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
% Translate to radial coordinates
[rows, ~] = size(TC.Position);
PolarPosition = zeros(rows, 2);
PolarVelocity = zeros(rows, 2);
Rin = 1.5;
Rout = 2.0;
win = 0.1;
wout = -0.075;

PolarPosition(:, 1) = sqrt(TC.Position(:, 1) .^ 2 + TC.Position(:, 2) .^ 2);
PolarPosition(:, 2) = atan2(TC.Position(:, 2), TC.Position(:, 1));

PolarVelocity(:, 1) = TC.Velocity(:, 1) .* TC.Position(:, 1) ./ PolarPosition(:, 1) + TC.Velocity(:, 2) .* TC.Position(:, 2) ./ PolarPosition(:, 1);
PolarVelocity(:, 2) = TC.Velocity(:, 2) .* TC.Position(:, 1) ./ PolarPosition(:, 1) - TC.Velocity(:, 1) .* TC.Position(:, 2) ./ PolarPosition(:, 1);

% Remove boundary particles
PolarPosition = PolarPosition(TC.Type == 1, :);
PolarVelocity = PolarVelocity(TC.Type == 1, :);

[ur, uth] = TC_Analytical(PolarPosition(:, 1), [Rin Rout win wout]);

figure;
plot(PolarPosition(:, 1), uth, 'o'); hold on;
plot(PolarPosition(:, 1), PolarVelocity(:, 2), '.');

L2_e = TC_Error(PolarPosition, PolarVelocity, ur, uth);

% fig1 = figure; hold on;
% TC.PlotParticles(fig1)

function [ur, uth] = TC_Analytical(r, params)
    Rin = params(1);
    Rout = params(2);
    win = params(3);
    wout = params(4);
    a =- ((Rout ^ 2 * Rin ^ 2) / (Rout ^ 2 - Rin ^ 2)) * (wout - win);
    b = (wout * Rout ^ 2 - win * Rin ^ 2) / (Rout ^ 2 - Rin ^ 2);
    ur = zeros(length(r), 1);
    uth = a ./ r + b * r;

end

function [L2] = TC_Error(PolarPosition, PolarVelocity, ur, uth)
    NormSPH = sqrt(PolarVelocity(:, 1) .^ 2 + PolarVelocity(:, 2) .^ 2);
    NormAnalytical = sqrt(ur .^ 2 + uth .^ 2);

    L2 = sqrt((1 / length(PolarPosition(:, 1))) * sum ((NormAnalytical - NormSPH) .^ 2));
end
