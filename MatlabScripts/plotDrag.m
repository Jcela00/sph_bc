clear all; close all; clc;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
% set(groot, 'defaultFigurePosition', [1440 0 600 600]);

% Load data

% accept any filename ended in DragLift.csv from parent folder ../
filename = dir('../*DragLift.csv').name;
filename = strcat('../', filename);
data = csvread(filename, 1, 0);

t = data(:, 1);
vx = data(:, 2);
vx2 = data(:, 3);
drag = data(:, 4);
drag2 = data(:, 5);
lift = data(:, 6);

t0 = 0.02/1.2e-4;
t = t / t0;
figure;
subplot(2, 1, 1);
hold on;
plot(t, drag);
plot(t, drag2);
yline(106.6, '--k')
axis([-inf inf 70 130]);
xlabel('Time [s]'); ylabel('$C_D$'); legend('Drag method 1', 'Drag method 2', 'Reference Value 106.6');

subplot(2, 1, 2); hold on;
plot(t, lift);
% axis([-inf inf 70 130]);
xlabel('Time [s]'); ylabel('$C_L$'); legend('Lift');

figure;
plot(t, vx); hold on;
plot(t, vx2);
yline(1.2e-4, '--k');
xlabel('Time [s]'); ylabel('$<u_x>$ [m/s]');
legend('$<u_x>$ method 1', '$<u_x>$ method 2', 'Reference 1.2e-4');
