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

rho0 = 1;
dim = 2;
Hconst = 1;

% filenames = {'../CSV_Data/Cavity/Cavity_new_summation_Re100_200_200_1rf_1prc/file',
%              '../CSV_Data/Cavity/Cavity_new_summation_Re1000_200_200_1rf_1prc/file',
%              '../CSV_Data/Cavity/Cavity_new_summation_Re10000_200_200_1rf_1prc/file'};

% names = {'Re100_200x200',
%          'Re1000_200x200'
%          'Re10000_200x200'};

% nfiles = [300 500 5000];

% figure; hold on;

% for k = 1:length(filenames)

%     [E stepvec] = CheckKineticEnergy(filenames{k}, 1, nfiles(k), rho0, dim, Hconst);
%     plot(stepvec, E, 'DisplayName', names{k});
% end

% save('CavityEnergy.mat', 'E');
% legend('Location', 'best');
% xlabel('$t$'); ylabel('$E$');
% exportgraphics(gcf, ['LatexFigures/CavityEnergy.pdf'], 'ContentType', 'vector', 'Resolution', 300);

lnwdth = 0.5;
mrksz = 2;

Npoints = 50;
xrange = [0 1];
yrange = [0 1];
streamDensity = 8;

filenames = {'../CSV_Data/Cavity/Cavity_new_summation_Re100_200_200_1rf_1prc/file',
             '../CSV_Data/Cavity/Cavity_new_summation_Re1000_200_200_1rf_1prc/file',
             '../CSV_Data/Cavity/Cavity_new_summation_Re10000_200_200_1rf_1prc/file', };

names = {'Re100_200x200',
         'Re1000_200x200',
         'Re10000_200x200'};

nfiles = [300 500 4500];

level1000 = [-3 -2 -1 -0.5 0 0.5 1 2 3 4 5];
level10000 = [-3 -2 0 1 2];

levelsVorticity = {level1000, level1000, level10000};
levels = [-1e-10 -1e-7 -1e-5 -1e-4 -0.01 -0.03 -0.05 -0.07 -0.09 -0.1 -0.11 -0.115 -0.1175 1e-8 1e-7 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];
levelsStreamfunction = {levels, levels, levels};

for k = 1:length(filenames)
    filename = filenames{k};
    Cavity = ParticleData(filename, nfiles(k), names{k}, rho0, dim, Hconst);

    fprintf('\nPlotting streamlines for %s\n', names{k});
    % fig1 = figure; hold on;
    fig2 = figure; hold on;
    fig3 = figure; hold on;
    fig4 = figure; hold on;
    Cavity.PlotStreamlines(fig2, fig3, fig4, lnwdth, xrange, yrange, Npoints, streamDensity, levelsStreamfunction{k}, levelsVorticity{k});

    % exportgraphics(fig1, ['LatexFigures/CavityStreamlines' num2str(Npoints) names{k} '.pdf'], 'ContentType', 'vector', 'Resolution', 300);
    exportgraphics(fig4, ['LatexFigures/CavityStreamfunction' num2str(Npoints) names{k} '.pdf'], 'ContentType', 'vector', 'Resolution', 300);
    exportgraphics(fig2, ['LatexFigures/CavityVorticityContour' num2str(Npoints) names{k} '.pdf'], 'ContentType', 'vector', 'Resolution', 300);
    exportgraphics(fig3, ['LatexFigures/CavityVelContour' num2str(Npoints) names{k} '.pdf'], 'ContentType', 'vector', 'Resolution', 300);

    close all;
end

function [E, stepvec] = CheckKineticEnergy(filename, nfiles_start, nfiles_end, rho0, dim, Hconst)

    E = [];
    step = 25;
    stepvec = [];

    for k = nfiles_start:step:nfiles_end

        if (k > 500)
            step = 100;
        end

        stepvec = [stepvec k];
        Cavity = ParticleData(filename, k, ['lolo'], rho0, dim, Hconst);
        vel = Cavity.Velocity(Cavity.FluidIndexes, 1:2);
        Etmp = vel(:, 1) .^ 2 + vel(:, 2) .^ 2;
        E = [E sum(Etmp) * Cavity.mass / 2];

    end

end
