clear all; close all; clc;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1.0);
set(groot, 'defaultLineMarkerSize', 2);
set(groot, 'defaultFigureUnits', 'centimeters');
set(groot, 'defaultFigurePosition', [100, 100, 16.0, 10.0]); %double column

% Reference time
t0 = 0.02/1.2e-4;

% Load data
dataset_sum_48 = ReadDragLift('../CSV_Data/Array/ArrayDrag/CylinderArray_new_summation_72_48_1rf_1prc_DragLift.csv', t0);
dataset_sum_96 = ReadDragLift('../CSV_Data/Array/ArrayDrag/CylinderArray_new_summation_144_96_1rf_1prc_DragLift.csv', t0);
dataset_sum_192 = ReadDragLift('../CSV_Data/Array/ArrayDrag/CylinderArray_new_summation_288_192_1rf_1prc_DragLift.csv', t0);

% Main plot
figure; hold on;
yline(106.76, '--k', 'DisplayName', 'Liu et al. (1998)');
plot(dataset_sum_48{1}, dataset_sum_48{3}, 'b-', 'DisplayName', '$N_y=48$');
plot(dataset_sum_96{1}, dataset_sum_96{3}, 'r-', 'DisplayName', '$N_y=96$');
plot(dataset_sum_192{1}, dataset_sum_192{3}, 'k-', 'DisplayName', '$N_y=192$');
xlabel('Time [s]');
ylabel('$C_D$');
axis([0 5 100 120]);
legend('Location', 'NorthWest', 'box', 'off');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);

% Define zoomed-in region
x_zoom = [2, 5];
y_zoom = [104, 108];

avg_48 = mean(dataset_sum_48{3}(dataset_sum_48{1} >= x_zoom(1) & dataset_sum_48{1} <= x_zoom(2)));
avg_96 = mean(dataset_sum_96{3}(dataset_sum_96{1} >= x_zoom(1) & dataset_sum_96{1} <= x_zoom(2)));
avg_192 = mean(dataset_sum_192{3}(dataset_sum_192{1} >= x_zoom(1) & dataset_sum_192{1} <= x_zoom(2)));
fprintf('Average drag coefficient for Ny=48: %.2f\n', avg_48);
fprintf('Average drag coefficient for Ny=96: %.2f\n', avg_96);
fprintf('Average drag coefficient for Ny=192: %.2f\n', avg_192);

% Add rectangle to indicate zoomed-in region
rectangle('Position', [x_zoom(1), y_zoom(1), diff(x_zoom), diff(y_zoom)], ...
    'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '--');

% Create the inset for the zoomed-in region
axesInset = axes('Position', [0.5, 0.5, 0.35, 0.35]); % Adjust inset position/size
hold on;
box on;
yline(106.76, '--k', 'DisplayName', 'Liu et al. (1998)');
plot(dataset_sum_48{1}, dataset_sum_48{3}, 'b-', 'DisplayName', '$N_y=48$');
plot(dataset_sum_96{1}, dataset_sum_96{3}, 'r-', 'DisplayName', '$N_y=96$');
plot(dataset_sum_192{1}, dataset_sum_192{3}, 'k-', 'DisplayName', '$N_y=192$');
xlim(x_zoom);
ylim(y_zoom);
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/ArrayDragCoefficient.pdf'], 'ContentType', 'vector', 'Resolution', 300);

% plot velocity
figure; hold on;
plot(dataset_sum_48{1}, dataset_sum_48{2}, 'b-', 'DisplayName', '$N_y=48$');
plot(dataset_sum_96{1}, dataset_sum_96{2}, 'r-', 'DisplayName', '$N_y=96$');
plot(dataset_sum_192{1}, dataset_sum_192{2}, 'k-', 'DisplayName', '$N_y=192$');

% Plot particles
rho0 = 1000; dim = 2; Hconst = 1;

mrksz = 2;
lnwdth = 1;

nfiles = [0 83 250 833];

for k = 1:length(nfiles)
    Array = ParticleData('../CSV_Data/Array/CylinderArray_new_summation_72_48_1rf_1prc/file', nfiles(k), ['lolo'], rho0, dim, Hconst);
    fig1 = figure; hold on;
    Array.PlotParticles(fig1, mrksz);
    xlabel('$x$'); ylabel('$y$');
    axis equal; axis tight;
    set(gca, 'FontSize', 11);
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
    exportgraphics(gcf, ['LatexFigures/ArrayParticles' num2str(nfiles(k)) '.pdf'], 'ContentType', 'vector', 'Resolution', 300);

end

function [DragDataset] = ReadDragLift(filename, t0)
    data = csvread(filename, 1, 0);
    t = data(:, 1) / t0;
    u = data(:, 2);
    drag = data(:, 3);
    lift = data(:, 4);
    DragDataset = {t, u, drag, lift};
end
