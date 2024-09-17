clear all; close all; clc;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 600 600]);

% Load data

% accept any filename ended in DragLift.csv from parent folder ../
filename = dir('../*DragLift.csv').name;
filename = strcat('../', filename);
data = csvread(filename, 1, 0);

t = data(:, 1);
v = data(:, 2);
drag = data(:, 3);
lift = data(:, 4);
vx = data(:, 5);
drag_alt = data(:, 6);
lift_alt = data(:, 7);

figure;
subplot(2, 1, 1); hold on;
plot(t, drag);
plot(t, drag_alt);
yline(106.6, '--k')
axis([-inf inf 70 130]);
xlabel('Time [s]'); ylabel('$C_D$'); legend('Drag', 'Drag Alt');

subplot(2, 1, 2); hold on;
plot(t, lift);
plot(t, lift_alt);
yline(106.6, '--k')
axis([-inf inf 70 130]);
xlabel('Time [s]'); ylabel('$C_L$'); legend('Lift', 'Lift Alt');

figure;
subplot(2, 1, 1);
plot(t, vx); hold on;
yline(1.2e-4, '--k');
xlabel('Time [s]'); ylabel('$<u_x>$ [m/s]');

subplot(2, 1, 2);
plot(t, v); hold on;
yline(1.2e-4, '--k');
xlabel('Time [s]'); ylabel('$<u>$ [m/s]');
