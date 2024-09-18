clear all; close all; clc;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 600 600]);

dx = 0.1;
kappa_max = 1 / dx;
kappa_max_2 = 1 / (2 * dx);
kappa_max_3 = 1 / (3 * dx);
kappa_grid = linspace(0, 2 * kappa_max, 1000);
vol = curvedVol(kappa_grid, dx);
vref = dx ^ 2;

figure; hold on;
plot(kappa_grid, vol(1, :), 'r-');
plot(kappa_grid, vol(2, :), 'g-');
plot(kappa_grid, vol(3, :), 'b-');
xline(kappa_max_3, '--', 'k=1/(3\Deltax)');
xline(kappa_max_2, '--', 'k=1/(2\Deltax)');
xline(kappa_max, '--', 'k=1/\Deltax');
yline(vref, '--', 'V_{flat}');
xlabel('$\kappa$');
ylabel('Volume');
legend('first', 'second', 'third');

r_grid = linspace(dx * 0.1, 30 * dx, 1000);
inv_r = 1 ./ r_grid;
vol = curvedVol(inv_r, dx);

figure; hold on;
plot(r_grid, vol(1, :), 'r-');
plot(r_grid, vol(2, :), 'g-');
plot(r_grid, vol(3, :), 'b-');
xline(dx, '--', '\Deltax');
xline(2 * dx, '--', '2\Deltax');
xline(3 * dx, '--', '3\Deltax');
yline(vref, '--', 'V_{flat}');
xlabel('$R_0$');
ylabel('Volume');

function [vol] = curvedVol(kappa, dx)
    ncol = size(kappa, 2);
    vol = zeros(3, ncol);

    for n = 1:3
        vol(n, :) = 0.5 * (2 * dx + dx * dx * kappa - 2 * n * dx * dx * kappa) * dx;
    end

    vol = max(0, vol);

end
