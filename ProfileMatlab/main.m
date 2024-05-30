clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);

% sph parameters
dp = 0.025;
rho0 = 1;
dim = 2;
filenum = 1500;
Nfluid = [40 1 60];
Nboundary = [3, 0, 0];

filenums = [1:10] * 1000;
fig4 = figure; hold on;
xfine = linspace(0, 1, 1000);
ufine = prof_a(xfine);
plot(xfine, ufine, 'k', 'DisplayName', 'Analytical');

for kk = 1:length(filenums)
    filenum = filenums(kk);
    partTest3 = ParticleData('CSV_Data/Vrel/file', filenum, ['Nu test' num2str(filenum)], dp, Nfluid, Nboundary, rho0, dim);

    obj = partTest3;
    %%% Plot ParaView like particles %%%
    % fig1 = figure;
    % obj.PlotParticles(fig1);
    % xline(0.75 - 2 * obj.H, 'r--', 'HandleVisibility', 'off');
    % xline(0.75 + 2 * obj.H, 'r--', 'HandleVisibility', 'off');

    %%% Plot Particle Scatter Plot %%%
    % fig3 = figure; hold on;
    % obj.ScatterParticlePlot(fig3)
    % xfine = linspace(0, 1, 1000);
    % ufine = prof_a(xfine);
    % plot(xfine, ufine, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Analytical');

    %%% Plot SPH interpolated profile %%%
    xsamples = 200; % sample points to evaluate SPH sum
    z_coord = 0.5; % channel z point to evaluate profile

    obj.PlotProfile(fig4, xsamples, z_coord, 0);
    legend;

    %%% Plot sum W for all particles & histogram %%%
    % checksum = obj.CheckKernel();
    % figure;
    % histogram(checksum, 50);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytical velocity profile

function u = prof_a(x)

    g = 0.1;
    nu = 0.0097325;
    L = 1;
    u = (g / (2 * nu)) .* x .* (L - x);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
