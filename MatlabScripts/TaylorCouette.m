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

nfile = 250;
test = ParticleData('../CSV_Data/TEST/file', nfile, ['80x1 TEST'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
tc_80_1 = ParticleData('../CSV_Data/TC_NEW/TaylorCouette_NEW_BC_Summation_80_80_4prc_1rf/file', nfile, ['80x1'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
tc_80_2 = ParticleData('../CSV_Data/TC_NEW/TaylorCouette_NEW_BC_Summation_80_80_4prc_2rf/file', nfile, ['80x2'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
tc_80_3 = ParticleData('../CSV_Data/TC_NEW/TaylorCouette_NEW_BC_Summation_80_80_4prc_3rf/file', nfile, ['80x3'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

tc_160_1 = ParticleData('../CSV_Data/TC_NEW/TaylorCouette_NEW_BC_Summation_160_160_4prc_1rf/file', nfile, ['160x1'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
tc_160_2 = ParticleData('../CSV_Data/TC_NEW/TaylorCouette_NEW_BC_Summation_160_160_4prc_2rf/file', nfile, ['160x2'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
tc_160_3 = ParticleData('../CSV_Data/TC_NEW/TaylorCouette_NEW_BC_Summation_160_160_4prc_3rf/file', nfile, ['160x3'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

tc_240_1 = ParticleData('../CSV_Data/TC_NEW/TaylorCouette_NEW_BC_Summation_240_240_4prc_1rf/file', nfile, ['240x1'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
tc_240_2 = ParticleData('../CSV_Data/TC_NEW/TaylorCouette_NEW_BC_Summation_240_240_4prc_2rf/file', nfile, ['240x2'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
tc_240_3 = ParticleData('../CSV_Data/TC_NEW/TaylorCouette_NEW_BC_Summation_240_240_4prc_3rf/file', nfile, ['240x3'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

tc_80_old = ParticleData('../CSV_Data/TC_NEW/TaylorCouette_OLD_BC_Summation_80_80_4prc/file', nfile, ['80 OLD'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
tc_160_old = ParticleData('../CSV_Data/TC_NEW/TaylorCouette_OLD_BC_Summation_160_160_4prc/file', nfile, ['160 OLD'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
tc_240_old = ParticleData('../CSV_Data/TC_NEW/TaylorCouette_OLD_BC_Summation_240_240_4prc/file', nfile, ['240 OLD'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

setTest = {test};
set80 = {tc_80_1, tc_80_2, tc_80_3};
set160 = {tc_160_1, tc_160_2, tc_160_3};
set240 = {tc_240_1, tc_240_2, tc_240_3};
set1 = {tc_80_1, tc_160_1, tc_240_1};
set2 = {tc_80_2, tc_160_2, tc_240_2};
set3 = {tc_80_3, tc_160_3, tc_240_3};
setold = {tc_80_old, tc_160_old, tc_240_old};

E_test = ComputeErrorAndPlotTC(setTest);
E_1 = ComputeErrorAndPlotTC(set1);
E_2 = ComputeErrorAndPlotTC(set2);
E_3 = ComputeErrorAndPlotTC(set3);
pause;
close all;
E_80 = ComputeErrorAndPlotTC(set80);
E_160 = ComputeErrorAndPlotTC(set160);
E_240 = ComputeErrorAndPlotTC(set240);
E_old = ComputeErrorAndPlotTC(setold);

figure;
res = [80 160 240];
loglog([res], E_1, 'o-'); hold on;
loglog(res, E_2, 'o-');
loglog(res, E_3, 'o-');
loglog([res], E_old, 'o-');
loglog([240], E_test, 'o-');
ylabel('$L_2$ error');
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
    errors = zeros(nres, 2);

    fig1 = figure;
    fig2 = figure;

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

    errors = errors(:, 2);
end
