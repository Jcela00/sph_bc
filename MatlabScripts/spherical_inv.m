clear all; close all; clc;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 1440 1440]);
dp = 0.025;
Hconst = 3;
Nfluid = [10 10 1];
Nboundary = [10 10 1];
rho0 = 1;
dim = 2;

triangle = ParticleData('../CSV_Data/tri_eq', 0, ['triangle'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

triangle = triangle.OnlyBoundaryFluid(0);

P_old = triangle.Position;
% remove walls at y= 0 and y=1
P_original = P_old(P_old(:, 2) > 0 & P_old(:, 2) < 1, :);

NN = 1;
H_big = 0.433013 + 2 * NN * dp;
D_big = H_big / cos(pi / 6); % cos(30)

rf = 1;
dx = dp / rf;

Nside = ceil(D_big / dx);

t_1 = [cos(pi / 6) sin(pi / 6) 0] * dx;
t_2 = [cos(pi / 6) -sin(pi / 6) 0] * dx;
t_3 = [0 1 0] * dx;

centres = zeros(3 * Nside - 3, 3);

Tip = [0.461325 - NN * dp 0.5 0];

UR = Tip + t_1 * (Nside - 1);
LR = Tip + t_2 * (Nside - 1);
idx = 1;

for k = 0:Nside - 1
    point = Tip + t_1 * k;
    centres(idx, :) = point;
    idx = idx + 1;
    % plot(point(1), point(2), 'ob');
end

for k = 1:Nside - 2
    point = UR - t_3 * k;
    centres(idx, :) = point;
    idx = idx + 1;
    % plot(point(1), point(2), '+g');
end

for k = 0:Nside - 2
    point = LR - t_2 * k;
    centres(idx, :) = point;
    idx = idx + 1;
    % plot(point(1), point(2), '*r');
end

ncentres = size(centres, 1);

for n = 1:ncentres
    centre = centres(n, :) + (1 - 2 * rand(1, 3)) * dp * 0.25;

    P_old = P_original;
    P_old = P_old - centre;

    figure; hold on; axis equal;
    plot(P_old(:, 1), P_old(:, 2), 'ob');
    plot(0, 0, 'xk');

    % only consider solid particles inside interaction radius of 3*dp
    Norm_mat = sqrt(sum(P_old .^ 2, 2));
    P_old = P_old(Norm_mat < 3 * dp, :);

    R = max(sqrt(sum((P_old) .^ 2, 2))); % radius is the maximum distance from the origin
    X = PlotCircle([0 0], R);
    X2 = PlotCircle([0 0], 3 * dp);
    plot(X(1, :), X(2, :), '-k');
    plot(X2(1, :), X2(2, :), '--k');

    P_hat = spherical_inversion(P_old, R);
    P_old = [0 0 0; P_old];
    P_hat = [0 0 0; P_hat];
    P_hat = P_hat(:, 1:2); % remove z component full of zeros
    k = convhull(P_hat);

    plot(P_hat(:, 1), P_hat(:, 2), 'or');
    plot(P_hat(k, 1), P_hat(k, 2), '-r');
    plot(P_old(k, 1), P_old(k, 2), '+b');
    title(['n = ', num2str(n)]);

    saveas(gcf, ['Figures/spherical_inv_', num2str(n), '.png']);
    close all;
end

function [P_new] = spherical_inversion(P_old, R)
    normP_matrix = sqrt(sum((P_old) .^ 2, 2));
    P_new = P_old + 2 * (R - normP_matrix) .* P_old ./ normP_matrix;
end

function [X] = PlotCircle(C, r)
    Nsamples = 100;
    theta = linspace(0, 2 * pi, Nsamples);
    X = [C(1) + r * cos(theta); C(2) + r * sin(theta)];
end
