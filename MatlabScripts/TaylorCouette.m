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

res = [120, 160, 240];
steps = 100:100;

E_1 = zeros(length(steps), length(res));

for k = 1:length(steps)
    nfile = steps(k);

    % tc_80_1 = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_NEW_BC_Summation_80_80_4prc_1rf/file', nfile, ['80x1'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    % tc_80_2 = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_NEW_BC_Summation_80_80_4prc_2rf/file', nfile, ['80x2'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    % tc_80_3 = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_NEW_BC_Summation_80_80_4prc_3rf/file', nfile, ['80x3'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    % tc_80_old = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_OLD_BC_Summation_80_80_4prc/file', nfile, ['80xold'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

    tc_120_1 = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_NEW_BC_Summation_120_120_4prc_1rf/file', nfile, ['120x1'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    tc_120_2 = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_NEW_BC_Summation_120_120_4prc_2rf/file', nfile, ['120x2'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    tc_120_3 = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_NEW_BC_Summation_120_120_4prc_3rf/file', nfile, ['120x3'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    tc_120_old = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_OLD_BC_Summation_120_120_4prc/file', nfile, ['120xold'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

    tc_160_1 = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_NEW_BC_Summation_160_160_4prc_1rf/file', nfile, ['160x1'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    tc_160_2 = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_NEW_BC_Summation_160_160_4prc_2rf/file', nfile, ['160x2'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    tc_160_3 = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_NEW_BC_Summation_160_160_4prc_3rf/file', nfile, ['160x3'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    tc_160_old = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_OLD_BC_Summation_160_160_4prc/file', nfile, ['160xold'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

    tc_240_1 = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_NEW_BC_Summation_240_240_4prc_1rf/file', nfile, ['240x1'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    tc_240_2 = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_NEW_BC_Summation_240_240_4prc_2rf/file', nfile, ['240x2'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    tc_240_3 = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_NEW_BC_Summation_240_240_4prc_3rf/file', nfile, ['240x3'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    tc_240_old = ParticleData('../CSV_Data/TC_NEWEST/TaylorCouette_OLD_BC_Summation_240_240_4prc/file', nfile, ['240xold'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

    % set_80 = {tc_80_1, tc_80_2, tc_80_3, tc_80_old};
    set_120 = {tc_120_1, tc_120_2, tc_120_3, tc_120_old};
    set_160 = {tc_160_1, tc_160_2, tc_160_3, tc_160_old};
    set_240 = {tc_240_1, tc_240_2, tc_240_3, tc_240_old};

    % E_80 = ComputeErrorAndPlotTC(set_80, 0);
    E_120 = ComputeErrorAndPlotTC(set_120, 0);
    E_160 = ComputeErrorAndPlotTC(set_160, 0);
    E_240 = ComputeErrorAndPlotTC(set_240, 0);

    E_1(k, :) = [E_120(1), E_160(1), E_240(1)];
    E_2(k, :) = [E_120(2), E_160(2), E_240(2)];
    E_3(k, :) = [E_120(3), E_160(3), E_240(3)];
    E_old(k, :) = [E_120(4), E_160(4), E_240(4)];

    k
end

figure; hold on;
semilogy(steps, E_1(:, 1), 'o-');
semilogy(steps, E_1(:, 2), 'o-');
semilogy(steps, E_1(:, 3), 'o-');
% semilogy(steps, E_1(:, 4), 'o-');
title('x1'); legend('80x80', '120x120', '160x160', '240x240');
figure; hold on;
semilogy(steps, E_old(:, 1), 'o-');
semilogy(steps, E_old(:, 2), 'o-');
semilogy(steps, E_old(:, 3), 'o-');
% semilogy(steps, E_old(:, 4), 'o-');
title('old'); legend('80x80', '120x120', '160x160', '240x240');

figure; hold on;
semilogy(steps, E_1(:, 1), 'o-');
semilogy(steps, E_2(:, 1), 'o-');
semilogy(steps, E_3(:, 1), 'o-');
title('80x80'); legend('x1', 'x2', 'x3');

E_1 = mean(E_1, 1);
E_2 = mean(E_2, 1);
E_3 = mean(E_3, 1);
E_old = mean(E_old, 1);

figure; hold on;
loglog(res, E_1, 'o-');
loglog(res, E_2, 'o-');
loglog(res, E_3, 'o-');
loglog(res, E_old, 'o-');
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

function [errors] = ComputeErrorAndPlotTC(datasets, plotsetting)

    nres = length(datasets);
    errors = zeros(nres, 2);

    Rin = 1.5; Rout = 2.0;
    win = 0.1; wout = -0.075;

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
        PolarVelocityT = zeros(rows, 2);

        % get position and velocity from dataset
        pos = TC.Position;
        vel = TC.Velocity;
        velT = TC.VelocityTransport;

        % Translate to polar coordinates 1 r, 2 theta
        PolarPosition(:, 1) = sqrt(pos(:, 1) .^ 2 + pos(:, 2) .^ 2);
        PolarPosition(:, 2) = atan2(pos(:, 2), pos(:, 1));
        PolarVelocity(:, 1) = vel(:, 1) .* pos(:, 1) ./ PolarPosition(:, 1) + vel(:, 2) .* pos(:, 2) ./ PolarPosition(:, 1);
        PolarVelocity(:, 2) = vel(:, 2) .* pos(:, 1) ./ PolarPosition(:, 1) - vel(:, 1) .* pos(:, 2) ./ PolarPosition(:, 1);
        PolarVelocityT(:, 1) = velT(:, 1) .* pos(:, 1) ./ PolarPosition(:, 1) + velT(:, 2) .* pos(:, 2) ./ PolarPosition(:, 1);
        PolarVelocityT(:, 2) = velT(:, 2) .* pos(:, 1) ./ PolarPosition(:, 1) - velT(:, 1) .* pos(:, 2) ./ PolarPosition(:, 1);

        % Remove boundary particles
        PolarPosition = PolarPosition(TC.Type == 1, :);
        PolarVelocity = PolarVelocity(TC.Type == 1, :);
        PolarVelocityT = PolarVelocityT(TC.Type == 1, :);

        [ur, uth] = TC_Analytical(PolarPosition(:, 1), [Rin Rout win wout]);

        errors(k, 1) = TC_Error(PolarPosition, PolarVelocity, ur, uth); % error in velocity
        errors(k, 2) = TC_Error(PolarPosition, PolarVelocityT, ur, uth); % error in transport velocity

        % fine grid for analytical solution plot
        rfine = linspace(Rin, Rout, 1000);
        [urfine, uthfine] = TC_Analytical(rfine, [Rin Rout win wout]);

        if (plotsetting == 1)
            figure(fig2);
            subplot(2, nres, k); hold on;
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

            % figure(fig1);
            % subplot(2, nres, k);
            % hold on;
            % plot(PolarPosition(:, 1), abs(PolarVelocity(:, 2) - uth), '.'); ylabel('$E_\theta$');
            % xlabel('$r$');
            % axis([1.5 2 0 0.1]);
            % title([TC.PlotName]);
            % subplot(2, nres, k + nres);
            % hold on;
            % plot(PolarPosition(:, 1), abs(PolarVelocity(:, 1) - ur), '.'); ylabel('$E_r$');
            % xlabel('$r$');
            % axis([1.5 2 0 0.1]);
            % title([" $L_2$ error: " + num2str(errors(k, 1)) "$L_2$ error (transport): " + num2str(errors(k, 2))])
        end

    end

    errors
    errors = errors(:, 1);
end
