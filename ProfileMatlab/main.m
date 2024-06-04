clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
% make default figures larger
% make figures appear in main monitor
set(groot, 'defaultFigurePosition', [1440 0 1440 1440]);

% sph parameters
dp = 0.025;
rho0 = 1;
dim = 2;
filenum_end = 202;
Nfluid = [40 1 60];
Nboundary = [3, 0, 0];
Hconst = sqrt(3);
% Analytical profile
xfine = linspace(0, 1, 1000);
ufine = prof_a_couette(xfine);
xsamples = 100; % sample points to evaluate SPH sum
z_coord = 0.5; % channel z point to evaluate profile

% LOAD DATA
Couette = ParticleData('CSV_Data/Couette/file', filenum_end, ['Couette'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
obj = Couette;

% PLOT OVER TIME
timesteps = [[1:2:20] 40:20:filenum_end];
fig1 = figure; hold on;
obj.PlotProfileOverTime(fig1, xsamples, z_coord, 0, timesteps);
plot(xfine, ufine, 'k', 'DisplayName', 'Analytical');
legend;

% PLOT JUST FINAL TIME
fig2 = figure; hold on;
obj.PlotProfile(fig2, xsamples, z_coord, 0);
plot(xfine, ufine, 'k', 'DisplayName', 'Analytical');
axis equal;

% % PARTICLE SCATTER PLOT
% fig3 = figure; hold on;
% obj.ScatterParticlePlot(fig2);
% plot(xfine, ufine, 'k', 'DisplayName', 'Analytical');
% legend;

% % PLOT LIKE PARAVIEW
% fig4 = figure;
% obj.PlotParticles(fig3);
% % CHECK KERNEL
% checksum = obj.CheckKernel();
% figure;
% histogram(checksum, 50);

% % % single filenum
% filenum = 20000;
% partTest3 = ParticleData(filename, filenum, ['Nu test' num2str(filenum)], dp, Nfluid, Nboundary, rho0, dim, Hconst);
% obj = partTest3;
% %% Plot ParaView like particles %%%
% % fig2 = figure;
% % obj.PlotParticles(fig2);
% % xline(0.75 - 2 * obj.H, 'r--', 'HandleVisibility', 'off');
% % xline(0.75 + 2 * obj.H, 'r--', 'HandleVisibility', 'off');

% %%% Plot Particle Scatter Plot %%%
% fig3 = figure; hold on;
% obj.ScatterParticlePlot(fig3)
% xfine = linspace(0, 1, 1000);
% ufine = prof_a(xfine);
% plot(xfine, ufine, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
% legend;

% %%% Plot SPH interpolated profile %%%
% % fig4 = figure; hold on;
% % xsamples = 100; % sample points to evaluate SPH sum
% % z_coord = 0.5; % channel z point to evaluate profile

% % obj.PlotProfile(fig4, xsamples, z_coord, 0);
% % plot(xfine, ufine, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
% % legend;

% % %%% Plot sum W for all particles & histogram %%%
% % checksum = obj.CheckKernel();
% % figure;
% % histogram(checksum, 50);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytical velocity profile

function u = prof_a_pouiseuille(x)

    g = 0.1;
    nu = 0.01;
    L = 1;
    u = (g / (2 * nu)) .* x .* (L - x);

end

function u = prof_a_couette(x)

    L = 1;
    utop = 1.25;

    u = utop * x / L;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
