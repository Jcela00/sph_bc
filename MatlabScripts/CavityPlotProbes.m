clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 600 600]);

probe0_40 = csvread('../CSV_Data/CAVITY/probes_0_Cavity_NEW_BC_Summation_40_40_4prc/file_500.csv', 1, 0);

vx_40 = probe0_40(:, 1);
pos1_40 = probe0_40(:, 3:5);
[pos1_40, idx] = sortrows(pos1_40);
vx_40 = vx_40(idx);

probe1_40 = csvread('../CSV_Data/CAVITY/probes_1_Cavity_NEW_BC_Summation_40_40_4prc/file_500.csv', 1, 0);

vy_40 = probe1_40(:, 1);
pos2_40 = probe1_40(:, 3:5);
[pos2_40, idx] = sortrows(pos2_40);
vy_40 = vy_40(idx);

probe0_80 = csvread('../CSV_Data/CAVITY/probes_0_Cavity_NEW_BC_Summation_80_80_4prc/file_500.csv', 1, 0);

vx_80 = probe0_80(:, 1);
pos1_80 = probe0_80(:, 3:5);
[pos1_80, idx] = sortrows(pos1_80);
vx_80 = vx_80(idx);

probe1_80 = csvread('../CSV_Data/CAVITY/probes_1_Cavity_NEW_BC_Summation_80_80_4prc/file_500.csv', 1, 0);

vy_80 = probe1_80(:, 1);
pos2_80 = probe1_80(:, 3:5);
[pos2_80, idx] = sortrows(pos2_80);
vy_80 = vy_80(idx);

% ghia data
vxref = csvread('../CSV_Data/CAVITY/Ghia/vx.csv', 1, 0);
vyref = csvread('../CSV_Data/CAVITY/Ghia/vy.csv', 1, 0);

t = tiledlayout(1, 1);
ax1 = axes(t);
h1 = plot(ax1, vyref(:, 1), vyref(:, 2), 'sr', 'DisplayName', 'Ghia, $V_y$, 257x257');
hold on;
h2 = plot(ax1, pos2_40(:, 1), vy_40, '--k', 'DisplayName', 'sph 40x40');
h3 = plot(ax1, pos2_80(:, 1), vy_80, '-k', 'DisplayName', 'sph 80x80');
xlabel('$x$'); ylabel('$V_y(x)$');
axis([0.0 1.10 -0.6 0.5]);

ax2 = axes(t);
h5 = plot(ax2, vxref(:, 2), vxref(:, 1), 'ob', 'DisplayName', 'Ghia, $V_x$, 257x257');
hold on;
h6 = plot(ax2, vx_40, pos1_40(:, 2), '--k');
h7 = plot(ax2, vx_80, pos1_80(:, 2), '-k');
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
l = legend([h1, h5, h2, h3], 'Location', 'SouthEast');
l.Box = 'off';
xlabel('$V_x(y)$'); ylabel('$y$');
axis([-0.4 1.0 0.0 1.0]);
