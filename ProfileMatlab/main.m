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
filenum_end = 1050;
Nfluid = [40 1 20];
Nboundary = [3, 0, 0];
Hconst = sqrt(3);

% Analytical profile
xfine = linspace(0, 1, 1000);
ufine_p = prof_a_pouiseuille(xfine, 1050);
ufine_c = prof_a_couette(xfine, 1050);
xsamples = 100; % sample points to evaluate SPH sum
z_coord = 0.1; % channel z point to evaluate profile

% LOAD DATA
% Poiseuille_NoP_ArtVisc = ParticleData('../CSV_Data/Poiseuille_NoP_ArtVisc/file', filenum_end, ['No Pressure Art Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
% Poiseuille_OldP_ArtVisc = ParticleData('../CSV_Data/Poiseuille_OldP_ArtVisc/file', filenum_end, ['Old Pressure Art Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
% Poiseuille_NewP_ArtVisc = ParticleData('../CSV_Data/Poiseuille_NewP_ArtVisc/file', filenum_end, ['New Pressure Art Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));

% Poiseuille_NoP_PhysVisc = ParticleData('../CSV_Data/Poiseuille_NoP_PhysVisc/file', filenum_end, ['No Pressure Phys Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
% Poiseuille_OldP_PhysVisc = ParticleData('../CSV_Data/Poiseuille_OldP_PhysVisc/file', filenum_end, ['Old Pressure Phys Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
% Poiseuille_NewP_PhysVisc = ParticleData('../CSV_Data/Poiseuille_NewP_PhysVisc/file', filenum_end, ['New Pressure Phys Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));

% Couette_NoP_ArtVisc = ParticleData('../CSV_Data/Couette_NoP_ArtVisc/file', filenum_end, ['No Pressure Art Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
% Couette_OldP_ArtVisc = ParticleData('../CSV_Data/Couette_OldP_ArtVisc/file', filenum_end, ['Old Pressure Art Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));
% Couette_NewP_ArtVisc = ParticleData('../CSV_Data/Couette_NewP_ArtVisc/file', filenum_end, ['New Pressure Art Visc'], dp, Nfluid, Nboundary, rho0, dim, sqrt(3));

Couette_PhysVisc_NoTensile_NewPressure_KickDrift_1proc = ParticleData('../CSV_Data/Couette_PhysVisc_NoTensile_NewPressure_KickDrift_1proc/file', filenum_end, ['No Tensile New Pressure Kick Drift 1proc'], dp, Nfluid, Nboundary, rho0, dim, 1);
obj1 = Couette_PhysVisc_NoTensile_NewPressure_KickDrift_1proc;

timesteps = [1:100:1050];
fig1 = figure; hold on;
obj1.PlotProfileOverTime(fig1, xsamples, z_coord, 0, timesteps);
plot(xfine, ufine_c, 'k', 'DisplayName', 'Analytical');
legend;

% for k = 1:length(timesteps)
%     tt = timesteps(k);
%     plot(xfine, prof_a_couette(xfine, tt), 'DisplayName', ['t = ' num2str(tt)]);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytical velocity profiles

function u = prof_a_pouiseuille(x, t)

    g = 0.1;
    nu = 0.01;
    L = 1;
    Nterms = 20;
    u = (g / (2 * nu)) .* x .* (L - x);

    for n = 0:Nterms
        u = u - ((4 * g * L ^ 2) / (nu * pi ^ 3 * (2 * n + 1) ^ 3)) .* sin(pi * x * (2 * n + 1) / L) .* exp(- ((2 * n + 1) ^ 2 * pi ^ 2 * nu * t) / (L ^ 2));
    end

end

function u = prof_a_couette(x, t)

    L = 1;
    V0 = 1.25;
    nu = 0.01;
    Nterms = 10;
    u = V0 * x / L;

    for n = 1:Nterms
        u = u + ((2 * V0 / (n * pi)) * (-1) ^ n) .* sin(n * pi * x / L) .* exp(-nu * n ^ 2 * pi ^ 2 * t / (L ^ 2));
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
