clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 600 600]);
dp = 0.05;
Nfluid = [80 80 1];
Nboundary = [0, 0, 0];
Hconst = 1;
dim = 2;
rho0 = 1;

nfile = 148;
TC80_1 = ParticleData('../CSV_Data/TaylorCouette_NEW_BC_Summation_80_80_4prc_1rf_/file', nfile, ['80x80 1x wall'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
TC80_2 = ParticleData('../CSV_Data/TaylorCouette_NEW_BC_Summation_80_80_4prc_2rf_/file', nfile, ['80x80 2x wall'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
TC80_3 = ParticleData('../CSV_Data/TaylorCouette_NEW_BC_Summation_80_80_4prc_3rf_/file', nfile, ['80x80 3x wall'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

TC120_1 = ParticleData('../CSV_Data/TaylorCouette_NEW_BC_Summation_120_120_4prc_1rf_/file', nfile, ['120x120 1x wall'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
TC120_2 = ParticleData('../CSV_Data/TaylorCouette_NEW_BC_Summation_120_120_4prc_2rf_/file', nfile, ['120x120 2x wall'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
TC120_3 = ParticleData('../CSV_Data/TaylorCouette_NEW_BC_Summation_120_120_4prc_3rf_/file', nfile, ['120x120 3x wall'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

TC240_1 = ParticleData('../CSV_Data/TaylorCouette_NEW_BC_Summation_240_240_4prc_1rf_/file', nfile, ['240x240 1x wall'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
TC240_2 = ParticleData('../CSV_Data/TaylorCouette_NEW_BC_Summation_240_240_4prc_2rf_/file', nfile, ['240x240 2x wall'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
TC240_3 = ParticleData('../CSV_Data/TaylorCouette_NEW_BC_Summation_240_240_4prc_3rf_/file', nfile, ['240x240 3x wall'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

TC80_OLD = ParticleData('../CSV_Data/TaylorCouette_OLD_BC_Summation_80_80_4prc_1rf_/file', nfile, ['80x80 old'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
TC120_OLD = ParticleData('../CSV_Data/TaylorCouette_OLD_BC_Summation_120_120_4prc_1rf_/file', nfile, ['120x120 old'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
TC240_OLD = ParticleData('../CSV_Data/TaylorCouette_OLD_BC_Summation_240_240_4prc_1rf/file', nfile, ['240x240 old'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

datasets = {TC80_1 TC80_2 TC80_3 TC120_1 TC120_2 TC120_3 TC240_1 TC240_2 TC240_3};
data_single = {TC80_1 TC120_1 TC240_1};
data_double = {TC80_2 TC120_2 TC240_2};
data_triple = {TC80_3 TC120_3 TC240_3};
datasets_old = {TC80_OLD TC120_OLD TC240_OLD};

data_80 = {TC80_1 TC80_2 TC80_3};
data_120 = {TC120_1 TC120_2 TC120_3};
data_240 = {TC240_1 TC240_2 TC240_3};

E_single = ComputeErrorAndPlotTC(data_single);
E_double = ComputeErrorAndPlotTC(data_double);
E_triple = ComputeErrorAndPlotTC(data_triple);

close all;

E_old = ComputeErrorAndPlotTC(datasets_old);
E_80 = ComputeErrorAndPlotTC(data_80);
E_120 = ComputeErrorAndPlotTC(data_120);
E_240 = ComputeErrorAndPlotTC(data_240);

figure;
res = [80 120 240];
loglog(res, E_single, 'o-'); hold on;
loglog(res, E_double, 'o-');
loglog(res, E_triple, 'o-');
loglog(res, E_old, 'o-'); ylabel('$L_2$ error');
xlabel('$N_{particles}$');
legend('New BC', 'New BC 2x refined wall', 'New BC 3x refined wall', 'Old BC', 'Location', 'best');

function [ur, uth] = TC_Analytical(r, params)
    Rin = params(1);
    Rout = params(2);
    win = params(3);
    wout = params(4);
    a =- ((Rout ^ 2 * Rin ^ 2) / (Rout ^ 2 - Rin ^ 2)) * (wout - win);
    b = (wout * Rout ^ 2 - win * Rin ^ 2) / (Rout ^ 2 - Rin ^ 2);
    ur = zeros(length(r), 1);
    uth = a ./ r + b * r;

end

function [L2] = TC_Error(PolarPosition, PolarVelocity, ur, uth)
    NormSPH = sqrt(PolarVelocity(:, 1) .^ 2 + PolarVelocity(:, 2) .^ 2);
    NormAnalytical = sqrt(ur .^ 2 + uth .^ 2);

    L2 = sqrt((1 / length(PolarPosition(:, 1))) * sum ((NormAnalytical - NormSPH) .^ 2));
end

function [L2] = R_Error(PolarPosition, PolarVelocity, ur, uth)
    NormSPH = sqrt(PolarVelocity(:, 1) .^ 2);
    NormAnalytical = sqrt(ur .^ 2);

    L2 = sqrt((1 / length(PolarPosition(:, 1))) * sum ((NormAnalytical - NormSPH) .^ 2));
end

function [L2] = TH_Error(PolarPosition, PolarVelocity, ur, uth)
    NormSPH = sqrt(PolarVelocity(:, 2) .^ 2);
    NormAnalytical = sqrt(uth .^ 2);

    L2 = sqrt((1 / length(PolarPosition(:, 1))) * sum ((NormAnalytical - NormSPH) .^ 2));
end

function [errors] = ComputeErrorAndPlotTC(datasets)

    nres = length(datasets);
    errors = zeros(nres, 1);

    figure;

    for k = 1:nres
        TC = datasets{k};

        % Translate to radial coordinates
        [rows, ~] = size(TC.Position);
        PolarPosition = zeros(rows, 2);
        PolarVelocity = zeros(rows, 2);
        Rin = 1.5;
        Rout = 2.0;
        win = 0.1;
        wout = -0.075;

        PolarPosition(:, 1) = sqrt(TC.Position(:, 1) .^ 2 + TC.Position(:, 2) .^ 2);
        PolarPosition(:, 2) = atan2(TC.Position(:, 2), TC.Position(:, 1));

        PolarVelocity(:, 1) = TC.Velocity(:, 1) .* TC.Position(:, 1) ./ PolarPosition(:, 1) + TC.Velocity(:, 2) .* TC.Position(:, 2) ./ PolarPosition(:, 1);
        PolarVelocity(:, 2) = TC.Velocity(:, 2) .* TC.Position(:, 1) ./ PolarPosition(:, 1) - TC.Velocity(:, 1) .* TC.Position(:, 2) ./ PolarPosition(:, 1);

        % Remove boundary particles
        PolarPosition = PolarPosition(TC.Type == 1, :);
        PolarVelocity = PolarVelocity(TC.Type == 1, :);

        [ur, uth] = TC_Analytical(PolarPosition(:, 1), [Rin Rout win wout]);

        rfine = linspace(Rin, Rout, 1000);
        [urfine, uthfine] = TC_Analytical(rfine, [Rin Rout win wout]);
        subplot(2, nres, k);
        hold on;
        plot(PolarPosition(:, 1), PolarVelocity(:, 2), '.'); ylabel('$u_\theta$');
        plot(rfine, uthfine);
        title(TC.PlotName);
        axis([1.5 2 -0.15 0.15]);
        subplot(2, nres, k + nres);
        hold on;
        plot(PolarPosition(:, 1), PolarVelocity(:, 1), '.'); ylabel('$u_r$');
        plot(rfine, urfine);
        axis([1.5 2 -0.1 0.1]);
        title(TC.PlotName);

        errors(k) = TC_Error(PolarPosition, PolarVelocity, ur, uth);
    end

end
