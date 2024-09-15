clear all; close all; clc;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 600 600]);

% Load data

data = csvread('../Drag_Lift.csv', 1, 0);

t = data(:, 1);
vx = data(:, 2);
vy = data(:, 3);
v = data(:, 4);
vtx = data(:, 5);
vty = data(:, 6);
vt = data(:, 7);
drag = data(:, 8);
lift = data(:, 9);

figure;
subplot(1, 2, 1);
plot(t, drag);
yline(106.6);
xlabel('Time [s]');
ylabel('Drag [N]'); axis([[-inf inf 50 200]]);
subplot(1, 2, 2);
plot(t, lift);
xlabel('Time [s]');
ylabel('Lift [N]');

figure;
plot(t, vx, 'DisplayName', 'vx'); hold on;
yline(1.2e-4, 'DisplayName', 'vx = 1.2e-4 m/s', 'LineWidth', 1.5);
xlabel('Time [s]'); legend;
