clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 1440 1440]);

% sph parameters
dp = 0.025;
rho0 = 1;
dim = 2;
filenum_end = 100;
Nfluid = [40 1 60];
Nboundary = [3, 0, 0];
Hconst = sqrt(3);

% Analytical profile
xfine = linspace(0, 1, 1000);
ufine_p = prof_a_pouiseuille(xfine);
ufine_c = prof_a_couette(xfine);
xsamples = 100; % sample points to evaluate SPH sum
z_coord = 0.5; % channel z point to evaluate profile

% LOAD DATA
Poiseuille_NoP_ArtVisc = ParticleData('../CSV_Data/Poiseuille_NoP_ArtVisc/file', filenum_end, ['No Pressure Art Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
Poiseuille_OldP_ArtVisc = ParticleData('../CSV_Data/Poiseuille_OldP_ArtVisc/file', filenum_end, ['Old Pressure Art Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
Poiseuille_NewP_ArtVisc = ParticleData('../CSV_Data/Poiseuille_NewP_ArtVisc/file', filenum_end, ['New Pressure Art Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));

Poiseuille_NoP_PhysVisc = ParticleData('../CSV_Data/Poiseuille_NoP_PhysVisc/file', filenum_end, ['No Pressure Phys Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
Poiseuille_OldP_PhysVisc = ParticleData('../CSV_Data/Poiseuille_OldP_PhysVisc/file', filenum_end, ['Old Pressure Phys Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
Poiseuille_NewP_PhysVisc = ParticleData('../CSV_Data/Poiseuille_NewP_PhysVisc/file', filenum_end, ['New Pressure Phys Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));

Couette_NoP_ArtVisc = ParticleData('../CSV_Data/Couette_NoP_ArtVisc/file', filenum_end, ['No Pressure Art Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
Couette_OldP_ArtVisc = ParticleData('../CSV_Data/Couette_OldP_ArtVisc/file', filenum_end, ['Old Pressure Art Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
Couette_NewP_ArtVisc = ParticleData('../CSV_Data/Couette_NewP_ArtVisc/file', filenum_end, ['New Pressure Art Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));

obj1 = Poiseuille_NoP_ArtVisc;
obj2 = Poiseuille_OldP_ArtVisc;
obj3 = Poiseuille_NewP_ArtVisc;

obj4 = Poiseuille_NoP_PhysVisc;
obj5 = Poiseuille_OldP_PhysVisc;
obj6 = Poiseuille_NewP_PhysVisc;

obj7 = Couette_NoP_ArtVisc;
obj8 = Couette_OldP_ArtVisc;
obj9 = Couette_NewP_ArtVisc;

% PLOT OVER TIME
% POUISEUILLE ART VISC
% timesteps = [[1:2:20] 40:20:filenum_end];
% fig1 = figure; hold on;
% obj1.PlotProfileOverTime(fig1, xsamples, z_coord, 0, timesteps);
% plot(xfine, ufine_p, 'k', 'DisplayName', 'Analytical');
% legend;

% timesteps = [[1:2:20] 40:20:filenum_end];
% fig1 = figure; hold on;
% obj2.PlotProfileOverTime(fig1, xsamples, z_coord, 0, timesteps);
% plot(xfine, ufine_p, 'k', 'DisplayName', 'Analytical');
% legend;

% timesteps = [[1:2:20] 40:20:filenum_end];
% fig1 = figure; hold on;
% obj3.PlotProfileOverTime(fig1, xsamples, z_coord, 0, timesteps);
% plot(xfine, ufine_p, 'k', 'DisplayName', 'Analytical');
% legend;

% % POUISEUILLE PHYS VISC
timesteps = [[1:2:20] 40:20:filenum_end];
fig1 = figure; hold on;
obj4.PlotProfileOverTime(fig1, xsamples, z_coord, 0, timesteps);
plot(xfine, ufine_p, 'k', 'DisplayName', 'Analytical');
legend;

% timesteps = [[1:2:20] 40:20:filenum_end];
% fig1 = figure; hold on;
% obj5.PlotProfileOverTime(fig1, xsamples, z_coord, 0, timesteps);
% plot(xfine, ufine_p, 'k', 'DisplayName', 'Analytical');
% legend;

% timesteps = [[1:2:20] 40:20:filenum_end];
% fig1 = figure; hold on;
% obj6.PlotProfileOverTime(fig1, xsamples, z_coord, 0, timesteps);
% plot(xfine, ufine_p, 'k', 'DisplayName', 'Analytical');
% legend;

% COUETTE ART VISC
%timesteps = [[1:2:20] 40:20:filenum_end];
%fig1 = figure; hold on;
%obj7.PlotProfileOverTime(fig1, xsamples, z_coord, 0, timesteps);
%plot(xfine, ufine_c, 'k', 'DisplayName', 'Analytical');
%legend;

%timesteps = [[1:2:20] 40:20:filenum_end];
%fig1 = figure; hold on;
%obj8.PlotProfileOverTime(fig1, xsamples, z_coord, 0, timesteps);
%plot(xfine, ufine_c, 'k', 'DisplayName', 'Analytical');
%legend;

%timesteps = [[1:2:20] 40:20:filenum_end];
%fig1 = figure; hold on;
%obj9.PlotProfileOverTime(fig1, xsamples, z_coord, 0, timesteps);
%plot(xfine, ufine_c, 'k', 'DisplayName', 'Analytical');
%legend;

% PLOT JUST FINAL TIME
% POUISEUILLE
fig2 = figure; hold on;
obj1.PlotProfile(fig2, xsamples, z_coord, 0);
obj2.PlotProfile(fig2, xsamples, z_coord, 0);
obj3.PlotProfile(fig2, xsamples, z_coord, 0);
obj4.PlotProfile(fig2, xsamples, z_coord, 0);
obj5.PlotProfile(fig2, xsamples, z_coord, 0);
obj6.PlotProfile(fig2, xsamples, z_coord, 0);
plot(xfine, ufine_p, 'k', 'DisplayName', 'Analytical');
axis equal; legend;

% % COUETTE

% fig3 = figure; hold on;
% obj7.PlotProfile(fig3, xsamples, z_coord, 0);
% obj8.PlotProfile(fig3, xsamples, z_coord, 0);
% obj9.PlotProfile(fig3, xsamples, z_coord, 0);
% plot(xfine, ufine_c, 'k', 'DisplayName', 'Analytical');
% axis equal; legend;

% % PARTICLE SCATTER PLOT
% fig3 = figure; hold on;
% obj.ScatterParticlePlot(fig2);
% plot(xfine, ufine, 'k', 'DisplayName', 'Analytical');
% legend;

% % PLOT LIKE PARAVIEW
% fig4 = figure;
% obj.PlotParticles(fig3);

% CHECK KERNEL
% checksum = obj1.CheckKernel();
% figure;
% histogram(checksum, 50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytical velocity profiles

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
