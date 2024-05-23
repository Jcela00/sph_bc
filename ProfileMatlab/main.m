clear all; clc; close all;

%
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);

% sph parameters
dp = 0.025;
rho0 = 1000;
dim = 2;
filenum = 1500;
Nfluid = [40 1 60];
Nboundary = [3, 0, 0];
%partFS = ParticleData('FS/free_slip', filenum, 'Free Slip', dp, Nfluid, Nboundary, rho0, dim);
%partNS = ParticleData('NS/no_slip', filenum, 'No Slip', dp, Nfluid, Nboundary, rho0, dim);
partCbar600 = ParticleData('Cbar600/cbar600', 1600, 'No Slip', dp, Nfluid, Nboundary, rho0, dim);
partCbar1200 = ParticleData('Cbar1200/cbar1200', 672, 'No Slip', dp, Nfluid, Nboundary, rho0, dim);
partCbar300 = ParticleData('Cbar300/cbar300', 2000, 'No Slip', dp, Nfluid, Nboundary, rho0, dim);
partTest = ParticleData('NuTest/nu_test', 2000, 'Nu test', dp, Nfluid, Nboundary, rho0, dim);

fig1 = figure;
partTest.PlotParticles(fig1);
xline(0.5 - 2 * partTest.H, 'r--', 'HandleVisibility', 'off');
xline(0.5 + 2 * partTest.H, 'r--', 'HandleVisibility', 'off');

fig3 = figure; hold on;
partTest.ScatterParticlePlot(fig3)
xfine = linspace(0, 1, 1000);
ufine = prof_a(xfine);
plot(xfine, ufine, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Analytical');

fig4 = figure; hold on;
xsamples = 100;
zsamples = 1;
partTest.PlotProfile(fig4, xsamples, zsamples, 0);
plot(xfine, ufine, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
legend;

checksum = partTest.CheckKernel();
figure;
histogram(checksum, 50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytical velocity profile

function u = prof_a(x)
    % u =g/(2*nu)*x*(h-x)
    % or umax= g*h^2/(8*nu) and then
    % u = (4umax/h)*x(1-x/h)

    % gravity: 9.81
    % Re: 184.752
    % L: 1
    % nu: 0.0814695
    % alpha: 0.1
    % umax: 15.0517
    % cbar: 150.517
    % B: 3.23646e+06

    g = 9.81;
    nu = 0.0433012701892219 * 150.517 * 0.1/8/2;
    L = 1;
    % w = (g / (2 * nu)) .* x .* (h - x);
    u = (g / (2 * nu)) .* x .* (L - x);
    % u = (4 * umax / h) * x .* (1 - x / h);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
