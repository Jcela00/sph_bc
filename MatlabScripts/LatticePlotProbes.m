clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 600 600]);
dp = 0.002;
Nfluid = [50 50 1];
Nboundary = [0, 1, 0];
Hconst = 1;
dim = 2;
rho0 = 1000;

latticeNew = ParticleData('../CSV_Data/LATTICE/CylinderLattice_NEW_BC_Summation_50_50_4prc/file', 1000, ['new 50x50'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
latticeNew.Velocity = latticeNew.Velocity / 1e-4;
latticeNew.VelMax = max(sqrt(sum(latticeNew.Velocity .^ 2, 2)));
latticeNew.VelMin = min(sqrt(sum(latticeNew.Velocity .^ 2, 2)));
levels = [0.2 0.4 0.6 0.8 1.0];
fig1 = figure; hold on;
latticeNew.PlotParticles(fig1);
PlotMorrisContours(fig1)
latticeNew.PlotContour(fig1, levels);
axis equal;

latticeNew100 = ParticleData('../CSV_Data/LATTICE/CylinderLattice_NEW_BC_Summation_100_100_4prc/file', 1000, ['new 100x100'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
latticeNew100.Velocity = latticeNew100.Velocity / 1e-4;
latticeNew100.VelMax = max(sqrt(sum(latticeNew100.Velocity .^ 2, 2)));
latticeNew100.VelMin = min(sqrt(sum(latticeNew100.Velocity .^ 2, 2)));
levels = [0.2 0.4 0.6 0.8 1.0];
fig2 = figure; hold on;
latticeNew100.PlotParticles(fig2);
PlotMorrisContours(fig2)
latticeNew100.PlotContour(fig2, levels);
axis equal;

latticeOld = ParticleData('../CSV_Data/LATTICE/CylinderLattice_OLD_BC_Summation_50_50_4prc/file', 1000, ['old 50x50'], dp, [100 100 1], Nboundary, rho0, dim, Hconst);
latticeOld.Velocity = latticeOld.Velocity / 1e-4;
latticeOld.VelMax = max(sqrt(sum(latticeOld.Velocity .^ 2, 2)));
latticeOld.VelMin = min(sqrt(sum(latticeOld.Velocity .^ 2, 2)));
levels = [0.2 0.4 0.6 0.8 1.0];
fig2 = figure; hold on;
latticeOld.PlotParticles(fig2);
PlotMorrisContours(fig2)
latticeOld.PlotContour(fig2, levels);
axis equal;

latticeOld_100 = ParticleData('../CSV_Data/LATTICE/CylinderLattice_OLD_BC_Summation_100_100_4prc/file', 1000, ['old 100x100'], dp, [100 100 1], Nboundary, rho0, dim, Hconst);
latticeOld_100.Velocity = latticeOld_100.Velocity / 1e-4;
latticeOld_100.VelMax = max(sqrt(sum(latticeOld_100.Velocity .^ 2, 2)));
latticeOld_100.VelMin = min(sqrt(sum(latticeOld_100.Velocity .^ 2, 2)));
levels = [0.2 0.4 0.6 0.8 1.0];
fig2 = figure; hold on;
latticeOld_100.PlotParticles(fig2);
PlotMorrisContours(fig2)
latticeOld_100.PlotContour(fig2, levels);
axis equal;

path1 = csvread('../CSV_Data/LATTICE/LatticeMorris/path1.csv', 1, 0);
path2 = csvread('../CSV_Data/LATTICE/LatticeMorris/path2.csv', 1, 0);

probe0 = csvread('../CSV_Data/LATTICE/probes_0_CylinderLattice_NEW_BC_Summation_50_50_4prc/file_1000.csv', 1, 0);
probe1 = csvread('../CSV_Data/LATTICE/probes_1_CylinderLattice_NEW_BC_Summation_50_50_4prc/file_1000.csv', 1, 0);

probe0_100 = csvread('../CSV_Data/LATTICE/probes_0_CylinderLattice_NEW_BC_Summation_100_100_4prc/file_1000.csv', 1, 0);
probe1_100 = csvread('../CSV_Data/LATTICE/probes_1_CylinderLattice_NEW_BC_Summation_100_100_4prc/file_1000.csv', 1, 0);

probe0_old = csvread('../CSV_Data/LATTICE/probes_0_CylinderLattice_OLD_BC_Summation_50_50_4prc/file_1000.csv', 1, 0);
probe1_old = csvread('../CSV_Data/LATTICE/probes_1_CylinderLattice_OLD_BC_Summation_50_50_4prc/file_1000.csv', 1, 0);

probe0_old_100 = csvread('../CSV_Data/LATTICE/probes_0_CylinderLattice_OLD_BC_Summation_100_100_4prc/file_1000.csv', 1, 0);
probe1_old_100 = csvread('../CSV_Data/LATTICE/probes_1_CylinderLattice_OLD_BC_Summation_100_100_4prc/file_1000.csv', 1, 0);

% map FEM coordinates from -0.5 0.5 to 0 0.01
path1(:, 1) = (path1(:, 1) + 0.5) / 10;
path2(:, 1) = (path2(:, 1) + 0.5) / 10;

% rescale by 1e-4
probe0(:, 1) = probe0(:, 1) / 1e-4;
probe1(:, 1) = probe1(:, 1) / 1e-4;

probe0_100(:, 1) = probe0_100(:, 1) / 1e-4;
probe1_100(:, 1) = probe1_100(:, 1) / 1e-4;

probe0_old(:, 1) = probe0_old(:, 1) / 1e-4;
probe1_old(:, 1) = probe1_old(:, 1) / 1e-4;

probe0_old_100(:, 1) = probe0_old_100(:, 1) / 1e-4;
probe1_old_100(:, 1) = probe1_old_100(:, 1) / 1e-4;

figure;
plot(path1(:, 1), path1(:, 2), 'k', 'DisplayName', 'Path 1 FEM'); hold on;
plot(path2(:, 1), path2(:, 2), 'k', 'DisplayName', 'Path 2 FEM');
plot(probe0(:, 4), probe0(:, 1), 'or', 'DisplayName', 'New SPH 50x50');
plot(probe1(:, 4), probe1(:, 1), 'or', 'HandleVisibility', 'off');
plot(probe0_100(:, 4), probe0_100(:, 1), 'sr', 'DisplayName', 'New SPH 100x100');
plot(probe1_100(:, 4), probe1_100(:, 1), 'sr', 'HandleVisibility', 'off');
plot(probe0_old(:, 4), probe0_old(:, 1), 'xb', 'DisplayName', 'Old SPH 50x50');
plot(probe1_old(:, 4), probe1_old(:, 1), 'xb', 'HandleVisibility', 'off');
plot(probe0_old_100(:, 4), probe0_old_100(:, 1), '+b', 'DisplayName', 'Old SPH 100x100');
plot(probe1_old_100(:, 4), probe1_old_100(:, 1), '+b', 'HandleVisibility', 'off');

xlabel('$x$');
ylabel('$y$');
legend('Location', 'best');

function PlotMorrisContours(fig)
    MorrisContours = csvread('../CSV_Data/LATTICE/LatticeMorris/contours.csv', 2, 0);
    figure(fig);

    for k = 1:10
        Contour = MorrisContours(:, 2 * k - 1:2 * k);
        Contour(~any(Contour, 2), :) = [];
        Contour = Contour / 10;
        plot(Contour(:, 1), Contour(:, 2), 'r', 'LineWidth', 4);
    end

end
