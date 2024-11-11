clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 6);
% This is for exportgraphics
set(groot, 'defaultFigureUnits', 'centimeters');
% set(groot, 'defaultFigurePosition', [0, 0, 8.5, 6.0]); %single column
set(groot, 'defaultFigurePosition', [0, 0, 16.0, 10.0]); %double column

rho0 = 1000;
dim = 2;
Hconst = 1;

lnwdth = 2;
mrksz = 2;

latticeNew = ParticleData('../CSV_Data/Lattice/CylinderLattice_new_summation_50_50_1rf_1prc/file', 100, ['new 50x50'], rho0, dim, Hconst);
Nfuid = length(latticeNew.FluidIndexes);
Nbound = length(latticeNew.ObstacleIndexes);
fprintf('50x50: number of fluid %d  number of boundary %d', Nfuid, Nbound);
latticeNew.Velocity = latticeNew.Velocity / 1e-4;
latticeNew.VelMax = max(sqrt(sum(latticeNew.Velocity .^ 2, 2)));
latticeNew.VelMin = min(sqrt(sum(latticeNew.Velocity .^ 2, 2)));
levels = [0.2 0.4 0.6 0.8 1.0];
fig1 = figure; hold on;
latticeNew.PlotParticles(fig1, mrksz);
PlotMorrisContours(fig1, lnwdth)
latticeNew.PlotVelocityContour(fig1, levels, lnwdth);
axis equal;
xlabel('$x$'); ylabel('$y$');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, 'LatexFigures/LatticeContour50.pdf', 'ContentType', 'vector', 'Resolution', 300);

latticeNew = ParticleData('../CSV_Data/Lattice/CylinderLattice_new_summation_100_100_1rf_1prc/file', 100, ['new 100x100'], rho0, dim, Hconst);
Nfluid = length(latticeNew.FluidIndexes);
Nbound = length(latticeNew.ObstacleIndexes);
fprintf('100x100: number of fluid %d  number of boundary %d', Nfluid, Nbound);
latticeNew.Velocity = latticeNew.Velocity / 1e-4;
latticeNew.VelMax = max(sqrt(sum(latticeNew.Velocity .^ 2, 2)));
latticeNew.VelMin = min(sqrt(sum(latticeNew.Velocity .^ 2, 2)));
levels = [0.2 0.4 0.6 0.8 1.0];
fig2 = figure; hold on;
latticeNew.PlotParticles(fig2, mrksz);
PlotMorrisContours(fig2, lnwdth)
latticeNew.PlotVelocityContour(fig2, levels, lnwdth);
axis equal;
xlabel('$x$'); ylabel('$y$');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, 'LatexFigures/LatticeContour100.pdf', 'ContentType', 'vector', 'Resolution', 300);

%%%%%%%% PLOT PROBES %%%%%%%%%
path1 = csvread('../CSV_Data/Lattice/Morris/path1.csv', 1, 0);
path2 = csvread('../CSV_Data/Lattice/Morris/path2.csv', 1, 0);

[vel0 y0] = ReadProbeData('../CSV_Data/Lattice/probes_0_CylinderLattice_new_summation_50_50_1rf_1prc/file_100.csv', 2);
[vel1 y1] = ReadProbeData('../CSV_Data/Lattice/probes_1_CylinderLattice_new_summation_50_50_1rf_1prc/file_100.csv', 2);

[vel0_100 y0_100] = ReadProbeData('../CSV_Data/Lattice/probes_0_CylinderLattice_new_summation_100_100_1rf_1prc/file_100.csv', 2);
[vel1_100 y1_100] = ReadProbeData('../CSV_Data/Lattice/probes_1_CylinderLattice_new_summation_100_100_1rf_1prc/file_100.csv', 2);

% map FEM coordinates from -0.5 0.5 to 0 0.01
path1(:, 1) = (path1(:, 1) + 0.5) / 10;
path2(:, 1) = (path2(:, 1) + 0.5) / 10;

% rescale by 1e-4
vel0(:, 1) = vel0(:, 1) / 1e-4;
vel1(:, 1) = vel1(:, 1) / 1e-4;

vel0_100(:, 1) = vel0_100(:, 1) / 1e-4;
vel1_100(:, 1) = vel1_100(:, 1) / 1e-4;

% for 50x50 take first half of the data
nelements = size(vel0, 1);
lastindex = floor(nelements / 2);
vel0 = vel0(1:lastindex); y0 = y0(1:lastindex);
vel1 = vel1(1:lastindex); y1 = y1(1:lastindex);

% for 100x100 take second half of the data
nelements = size(vel0_100, 1);
firstindex = floor(nelements / 2) +1;
vel0_100 = vel0_100(firstindex:end); y0_100 = y0_100(firstindex:end);
vel1_100 = vel1_100(firstindex:end); y1_100 = y1_100(firstindex:end);

figure;
plot(path1(:, 1), path1(:, 2), ':k', 'DisplayName', 'Path 1 FEM, Morris et al.'); hold on;
plot(path2(:, 1), path2(:, 2), 'k', 'DisplayName', 'Path 2 FEM, Morris et al.');
plot(y0, vel0, 'or', 'DisplayName', 'SPH 50x50');
plot(y1, vel1, 'or', 'HandleVisibility', 'off');
plot(y0_100, vel0_100, 'ob', 'DisplayName', 'SPH 100x100');
plot(y1_100, vel1_100, 'ob', 'HandleVisibility', 'off');

xlabel('$y$'); ylabel('$u_x$');
legend('Location', 'north', 'box', 'off');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, 'LatexFigures/LatticeProbes.pdf', 'ContentType', 'vector', 'Resolution', 300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotMorrisContours(fig, lnwdth)
    MorrisContours = csvread('../CSV_Data/Lattice/Morris/contours.csv', 2, 0);
    figure(fig);

    for k = 1:10
        Contour = MorrisContours(:, 2 * k - 1:2 * k);
        Contour(~any(Contour, 2), :) = [];
        Contour = Contour / 10;
        plot(Contour(:, 1), Contour(:, 2), 'r', 'LineWidth', lnwdth);
    end

end

function [velocity, position] = ReadProbeData(filename, position_component)
    % position_component: 1 -> x, 2 -> y, 3 -> z
    probeData = csvread(filename, 1, 0);

    velocity = probeData(:, 2);
    position = probeData(:, 4:6);
    position = position(:, position_component);
    [position, idx] = sortrows(position);
    velocity = velocity(idx);
end
