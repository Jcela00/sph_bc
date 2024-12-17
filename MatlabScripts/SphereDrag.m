clear all; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 2);
% This is for exportgraphics
set(groot, 'defaultFigureUnits', 'centimeters');
% set(groot, 'defaultFigurePosition', [0, 0, 8.5, 6.0]); %single column
set(groot, 'defaultFigurePosition', [100, 100, 16.0, 10.0]); %double column

U = 0.05;

D = 1;
R = D / 2;
A = pi * R * R;
params = [U A];
Re = U * D / 1;

dragData = ReadDragLiftVariableU('../CSV_Data/Sphere_new_summation_80_40_40_1rf_1prc/Sphere_new_summation_80_40_40_1rf_1prc_DragLift.csv', params);
stride = 100;
dragData{1} = dragData{1}(1:stride:end, :);
dragData{2} = dragData{2}(1:stride:end, :);
dragData{3} = dragData{3}(1:stride:end, :);
dragData{4} = dragData{4}(1:stride:end, :);

Re = dragData{2} * D / 1;

figure;
hold on;
plot(Re, 24 ./ Re, ':k', 'LineWidth', 3, 'DisplayName', 'Stokes Law');
plot(Re, dragData{3}, 'r', 'DisplayName', 'SPH');
legend('Location', 'Best');
axis([0 inf 0 4500]);
xlabel('$Re$');
ylabel('$C_d$');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, 'LatexFigures/SphereDrag.pdf', 'ContentType', 'vector', 'Resolution', 300);

Re_final = Re(end, :);
Stokes = 24 ./ Re_final;

figure; hold on;
yline(Stokes, 'DisplayName', 'Stokes');
plot(dragData{1}, dragData{3}, 'DisplayName', 'SPH');
legend('Location', 'Best');
axis([-inf inf 0 500]);
xlabel('$t$');
ylabel('$C_d$');

function [DragDataset] = ReadDragLift(filename, params)
    data = csvread(filename, 1, 0);
    t = data(:, 1);
    u = data(:, 2);
    drag = data(:, 3);
    lift = data(:, 4);

    U = params(1);
    A = params(2);
    drag = drag / (0.5 * U * U * A);
    lift = lift / (0.5 * U * U * A);

    DragDataset = {t, u, drag, lift};

end

function [DragDataset] = ReadDragLiftVariableU(filename, params)
    data = csvread(filename, 1, 0);
    t = data(:, 1);
    u = data(:, 2);
    drag = data(:, 3);
    lift = data(:, 4);

    U = params(1);
    A = params(2);
    drag = drag ./ (0.5 * u .* u * A);
    lift = lift ./ (0.5 * u .* u * A);

    DragDataset = {t, u, drag, lift};

end

function Cd = MorrisonEquation(Re)

    Cd = 24 ./ Re + (2.6 * Re / 5.0) ./ (1 + (Re / 5.0) .^ 1.52) + (0.411 * (Re / 263000) .^ -7.94) ./ (1 + (Re / 263000) .^ -8.00) + (0.25 * Re / 1e6) ./ (1 + Re / 1e6);
end
