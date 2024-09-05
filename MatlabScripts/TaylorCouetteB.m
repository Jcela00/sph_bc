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

nfile = 400;

steps = 1:400;
E1 = zeros(length(steps), 1);
E2 = zeros(length(steps), 1);
E3 = zeros(length(steps), 1);
E4 = zeros(length(steps), 1);
E5 = zeros(length(steps), 1);
E6 = zeros(length(steps), 1);

for k = 1:length(steps)
    nfile = steps(k);

    B1 = ParticleData('../CSV_Data/TC_Btest/TaylorCouette_NEW_BC_Summation_80_80_4prc_1rf_B1/file', nfile, ['B=1'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    B2 = ParticleData('../CSV_Data/TC_Btest/TaylorCouette_NEW_BC_Summation_80_80_4prc_1rf_B2/file', nfile, ['B=2'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    B3 = ParticleData('../CSV_Data/TC_Btest/TaylorCouette_NEW_BC_Summation_80_80_4prc_1rf_B3/file', nfile, ['B=3'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    B4 = ParticleData('../CSV_Data/TC_Btest/TaylorCouette_NEW_BC_Summation_80_80_4prc_1rf_B4/file', nfile, ['B=4'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    B5 = ParticleData('../CSV_Data/TC_Btest/TaylorCouette_NEW_BC_Summation_80_80_4prc_1rf_B5/file', nfile, ['B=5'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    B6 = ParticleData('../CSV_Data/TC_Btest/TaylorCouette_NEW_BC_Summation_80_80_4prc_1rf_B6/file', nfile, ['B=6'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

    set1 = {B1}; set2 = {B2}; set3 = {B3}; set4 = {B4}; set5 = {B5}; set6 = {B6};
    E1(k) = ComputeErrorAndPlotTC(set1, 0);
    E2(k) = ComputeErrorAndPlotTC(set2, 0);
    E3(k) = ComputeErrorAndPlotTC(set3, 0);
    E4(k) = ComputeErrorAndPlotTC(set4, 0);
    E5(k) = ComputeErrorAndPlotTC(set5, 0);
    E6(k) = ComputeErrorAndPlotTC(set6, 0);

    k
end

figure; hold on;
plot(steps, E1, '-');
plot(steps, E2, '-');
plot(steps, E3, '-');
plot(steps, E4, '-');
plot(steps, E5, '-');
plot(steps, E6, '-');
legend('B=1', 'B=2', 'B=3', 'B=4', 'B=5', 'B=6');

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

function [errors] = ComputeErrorAndPlotTC(datasets, plotsetting)

    nres = length(datasets);
    errors = zeros(nres, 2);

    if (plotsetting == 1)
        fig1 = figure;
        fig2 = figure;
    end

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

        vel = TC.Velocity;
        % vel = TC.VelocityTransport;

        PolarVelocity(:, 1) = vel(:, 1) .* TC.Position(:, 1) ./ PolarPosition(:, 1) + vel(:, 2) .* TC.Position(:, 2) ./ PolarPosition(:, 1);
        PolarVelocity(:, 2) = vel(:, 2) .* TC.Position(:, 1) ./ PolarPosition(:, 1) - vel(:, 1) .* TC.Position(:, 2) ./ PolarPosition(:, 1);

        % Remove boundary particles

        PolarPosition = PolarPosition(TC.Type == 1, :);
        PolarVelocity = PolarVelocity(TC.Type == 1, :);

        [ur, uth] = TC_Analytical(PolarPosition(:, 1), [Rin Rout win wout]);

        rfine = linspace(Rin, Rout, 1000);
        [urfine, uthfine] = TC_Analytical(rfine, [Rin Rout win wout]);

        errors(k, 1) = TC_Error(PolarPosition, PolarVelocity, ur, uth);

        %% SAME WITH TRANSPORT VELOCITY
        vel = TC.VelocityTransport;
        PolarPosition = zeros(rows, 2);
        PolarVelocityT = zeros(rows, 2);
        PolarPosition(:, 1) = sqrt(TC.Position(:, 1) .^ 2 + TC.Position(:, 2) .^ 2);
        PolarPosition(:, 2) = atan2(TC.Position(:, 2), TC.Position(:, 1));

        PolarVelocityT(:, 1) = vel(:, 1) .* TC.Position(:, 1) ./ PolarPosition(:, 1) + vel(:, 2) .* TC.Position(:, 2) ./ PolarPosition(:, 1);
        PolarVelocityT(:, 2) = vel(:, 2) .* TC.Position(:, 1) ./ PolarPosition(:, 1) - vel(:, 1) .* TC.Position(:, 2) ./ PolarPosition(:, 1);
        % Remove boundary particles
        PolarPosition = PolarPosition(TC.Type == 1, :);
        PolarVelocityT = PolarVelocityT(TC.Type == 1, :);

        errors(k, 2) = TC_Error(PolarPosition, PolarVelocityT, ur, uth);

        if (plotsetting == 1)

            figure(fig2);
            subplot(2, nres, k);
            hold on;
            plot(PolarPosition(:, 1), PolarVelocity(:, 2), '.'); ylabel('$u_\theta$');
            plot(rfine, uthfine);
            title([TC.PlotName]);
            axis([1.5 2 -0.15 0.15]);
            subplot(2, nres, k + nres);
            hold on;
            plot(PolarPosition(:, 1), PolarVelocity(:, 1), '.'); ylabel('$u_r$');
            plot(rfine, urfine);
            axis([1.5 2 -0.1 0.1]);
            title([" $L_2$ error: " + num2str(errors(k, 1)) "$L_2$ error (transport): " + num2str(errors(k, 2))])

            figure(fig1);
            subplot(2, nres, k);
            hold on;
            plot(PolarPosition(:, 1), abs(PolarVelocity(:, 2) - uth), '.'); ylabel('$E_\theta$');
            xlabel('$r$');
            axis([1.5 2 0 0.1]);
            title([TC.PlotName]);
            subplot(2, nres, k + nres);
            hold on;
            plot(PolarPosition(:, 1), abs(PolarVelocity(:, 1) - ur), '.'); ylabel('$E_r$');
            xlabel('$r$');
            axis([1.5 2 0 0.1]);
            title([" $L_2$ error: " + num2str(errors(k, 1)) "$L_2$ error (transport): " + num2str(errors(k, 2))])
        end

    end

    errors = errors(:, 2);
end
