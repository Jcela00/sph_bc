clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 1440 1440]);
Nfluid = [80 80 1];
Nboundary = [1, 1 1];
rho0 = 1;
dim = 2;
Hconst = 1;
dp = 0.05;

test = ParticleData('../CSV_Data/TC/TaylorCouette_new_summation_10_10_1rf_1prc/file', 300, ['this is a test'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

fig = figure; hold on;
test.PlotParticles(fig)
