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
filenum = 400;
Nfluid = [40 1 60];
Nboundary = [3, 0, 0];
Hconst = sqrt(3);

% LOAD DATA
NoPressureOldVisc = ParticleData('CSV_Data/NoPressureOldVisc/file', filenum, ['NoPressure'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
NoPresureNewVisc = ParticleData('CSV_Data/NoPressureNewVisc/file', filenum, ['NoPressure'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
obj = NoPresureNewVisc;
% PLOT OVER TIME
fig1 = figure; hold on;
xsamples = 100; % sample points to evaluate SPH sum
z_coord = 0.5; % channel z point to evaluate profile
obj.PlotProfileOverTime(fig1, xsamples, z_coord, 0, [10:50:400]);
xfine = linspace(0, 1, 1000);
ufine = prof_a(xfine);
plot(xfine, ufine, 'k', 'DisplayName', 'Analytical');
legend;

% PARTICLE SCATTER PLOT
fig2 = figure; hold on;
obj.ScatterParticlePlot(fig2);
plot(xfine, ufine, 'k', 'DisplayName', 'Analytical');
legend;

% PLOT LIKE PARAVIEW
fig3 = figure;
obj.PlotParticles(fig3);
% CHECK KERNEL
checksum = obj.CheckKernel();
figure;
histogram(checksum, 50);

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

function u = prof_a(x)

    g = 0.1;
    nu = 0.01;
    L = 1;
    u = (g / (2 * nu)) .* x .* (L - x);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
