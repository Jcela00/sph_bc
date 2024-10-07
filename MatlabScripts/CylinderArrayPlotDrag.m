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

data_72 = csvread("CSV_CylinderArray_drag/CylinderArray_new_summation_72_48_1rf_1prc_DragLift.csv", 1, 0);
data_144 = csvread("CSV_CylinderArray_drag/CylinderArray_new_summation_144_96_1rf_1prc_DragLift.csv", 1, 0);

datasets = {data_72, data_144};

fig1 = figure; hold on; yline(106.76, '--k')
xlabel('Time [s]'); ylabel('$C_D$'); axis([-inf inf 100 120]);

fig2 = figure; hold on;
xlabel('Time [s]'); ylabel('$C_L$'); axis([-inf inf -1 1]);

fig3 = figure; hold on; yline(1.2e-4, '--k');
xlabel('Time [s]'); ylabel('$<u_x>$ [m/s]'); axis([-inf inf 1.1e-4 1.3e-4]);

for k = 1:length(datasets)
    data = datasets{k};

    t = data(:, 1);
    vx = data(:, 2);
    drag = data(:, 3);
    lift = data(:, 4);

    t0 = 0.02/1.2e-4;
    t = t / t0;

    figure(fig1);
    plot(t, drag);

    figure(fig2);
    plot(t, lift);

    figure(fig3);
    plot(t, vx);

    steadyvx = mean(vx((ceil(0.5 * length(vx)):end)));
    text(t(end), vx(end), num2str(steadyvx), 'HorizontalAlignment', 'left');
end

figure(fig1);
legend('Liu et al.(1998)', '72x48', '144x96', 'Location', 'Best');

figure(fig2);
legend('72x48', '144x96', 'Location', 'Best');

figure(fig3);
legend('Reference', '72x48', '144x96', 'Location', 'Best');
