clear all; clc; close all;
global TypeBoundary TypeObstacle TypeFluid

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

res = [10, 20, 30, 40];
% res = 0.5 ./ res;
steps = [1:300];
dirname = 'TC';

E_1 = zeros(length(steps), length(res));
E_2 = zeros(length(steps), length(res));
% E_3 = zeros(length(steps), length(res));
% E_old = zeros(length(steps), length(res));

% Total number of iterations
nn = length(steps);
h = waitbar(0, 'Please wait...');

for k = 1:length(steps)
    nfile = steps(k);

    tc_10_1 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_10_10_1rf_1prc/file'], nfile, ['10x1'], rho0, dim, Hconst);
    tc_10_2 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_10_10_2rf_1prc/file'], nfile, ['10x2'], rho0, dim, Hconst);
    % tc_10_3 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_10_10_3rf_1prc/file'], nfile, ['10x3'], rho0, dim, Hconst);
    % tc_10_old = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_old_summation_10_10_1rf_1prc/file'], nfile, ['10 old'], rho0, dim, Hconst);

    % tc_15_1 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_15_15_1rf_1prc/file'], nfile, ['15x1'], rho0, dim, Hconst);
    % tc_15_2 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_15_15_2rf_1prc/file'], nfile, ['15x2'], rho0, dim, Hconst);
    % tc_15_3 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_15_15_3rf_1prc/file'], nfile, ['15x3'], rho0, dim, Hconst);
    % tc_15_old = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_old_summation_15_15_1rf_1prc/file'], nfile, ['15 old'], rho0, dim, Hconst);

    tc_20_1 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_20_20_1rf_1prc/file'], nfile, ['20x1'], rho0, dim, Hconst);
    tc_20_2 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_20_20_2rf_1prc/file'], nfile, ['20x2'], rho0, dim, Hconst);
    % tc_20_3 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_20_20_3rf_1prc/file'], nfile, ['20x3'], rho0, dim, Hconst);
    % tc_20_old = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_old_summation_20_20_1rf_1prc/file'], nfile, ['20 old'], rho0, dim, Hconst);

    tc_30_1 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_30_30_1rf_1prc/file'], nfile, ['30x1'], rho0, dim, Hconst);
    tc_30_2 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_30_30_2rf_1prc/file'], nfile, ['30x2'], rho0, dim, Hconst);
    % tc_30_3 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_30_30_3rf_1prc/file'], nfile, ['30x3'], rho0, dim, Hconst);
    % tc_30_old = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_old_summation_30_30_1rf_1prc/file'], nfile, ['30 old'], rho0, dim, Hconst);

    tc_40_1 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_40_40_1rf_1prc/file'], nfile, ['40x1'], rho0, dim, Hconst);
    tc_40_2 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_40_40_2rf_1prc/file'], nfile, ['40x2'], rho0, dim, Hconst);
    % tc_40_3 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_40_40_3rf_1prc/file'], nfile, ['40x3'], rho0, dim, Hconst);
    % tc_40_old = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_old_summation_40_40_1rf_1prc/file'], nfile, ['40 old'], rho0, dim, Hconst);

    set_10 = {tc_10_1, tc_10_2}; %, tc_10_3, tc_10_old};
    % set_15 = {tc_15_1, tc_15_2};%, tc_15_3, tc_15_old};
    set_20 = {tc_20_1, tc_20_2}; %, tc_20_3, tc_20_old};
    set_30 = {tc_30_1, tc_30_2}; % tc_30_3};
    set_40 = {tc_40_1, tc_40_2};

    plotsetting = 0;

    if k == length(steps)
        plotsetting = 1;
    end

    % Error of the 10, 20, 30 resolutions at time step nfile
    E_10 = ComputeErrorAndPlotTC(set_10, plotsetting);
    % E_15 = ComputeErrorAndPlotTC(set_15, 0);
    E_20 = ComputeErrorAndPlotTC(set_20, plotsetting);
    E_30 = ComputeErrorAndPlotTC(set_30, plotsetting);
    E_40 = ComputeErrorAndPlotTC(set_40, plotsetting);

    % Accumulate Error of the different x1 x2 x3 refinements (and old) for the different resolutions
    E_1(k, :) = [E_10(1), E_20(1), E_30(1), E_40(1)];
    E_2(k, :) = [E_10(2), E_20(2), E_30(2), E_40(2)];
    % E_3(k, :) = [E_10(3), E_20(3), E_30(3), E];
    % E_old(k, :) = [E_10(4), E_15(4), E_20(4), 0];

    disp(['Step: ' num2str(nfile)]);
    waitbar(k / nn, h, sprintf('Progress: %d percent ..', round((k / nn) * 100)));

end

close(h);

figure; hold on;
semilogy(steps, E_1(:, 1), 'o-');
semilogy(steps, E_1(:, 2), 'o-');
semilogy(steps, E_1(:, 3), 'o-');
semilogy(steps, E_1(:, 4), 'o-');
title('x1'); legend('10', '20', '30', '40');

figure; hold on;
semilogy(steps, E_2(:, 1), 'o-');
semilogy(steps, E_2(:, 2), 'o-');
semilogy(steps, E_2(:, 3), 'o-');
semilogy(steps, E_2(:, 4), 'o-');
title('x2'); legend('10', '20', '30', '40');

% figure; hold on;
% semilogy(steps, E_3(:, 1), 'o-');
% semilogy(steps, E_3(:, 2), 'o-');
% semilogy(steps, E_3(:, 3), 'o-');
% % semilogy(steps, E_3(:, 4), 'o-');
% title('x3'); legend('10', '20', '30');

% figure; hold on;
% semilogy(steps, E_old(:, 1), 'o-');
% semilogy(steps, E_old(:, 2), 'o-');
% semilogy(steps, E_old(:, 3), 'o-');
% % semilogy(steps, E_old(:, 4), 'o-');
% title('old'); legend('10', '20', '30');

figure; hold on;
semilogy(steps, E_1(:, 1), 'o-');
semilogy(steps, E_2(:, 1), 'o-');
% semilogy(steps, E_3(:, 1), 'o-');
title('10'); legend('x1', 'x2');
figure; hold on;

semilogy(steps, E_1(:, 2), 'o-');
semilogy(steps, E_2(:, 2), 'o-');
% semilogy(steps, E_3(:, 2), 'o-');
title('20'); legend('x1', 'x2');

figure; hold on;
semilogy(steps, E_1(:, 3), 'o-');
semilogy(steps, E_2(:, 3), 'o-');
% semilogy(steps, E_3(:, 3), 'o-');
title('30'); legend('x1', 'x2');

figure; hold on;
semilogy(steps, E_1(:, 4), 'o-');
semilogy(steps, E_2(:, 4), 'o-');
% semilogy(steps, E_3(:, 3), 'o-');
title('40'); legend('x1', 'x2');

averaging_window = [length(steps) - 1, length(steps)];
E_1_avg = mean(E_1(averaging_window(1):averaging_window(2), :), 1);
E_2_avg = mean(E_2(averaging_window(1):averaging_window(2), :), 1);
% E_3_avg = mean(E_3(averaging_window(1):averaging_window(2), :), 1);
% E_old_avg = mean(E_old(averaging_window(1):averaging_window(2), :), 1);

p1 = polyfit(log(res(1:4)), log(E_1_avg(1:4)), 1);
p2 = polyfit(log(res(1:4)), log(E_2_avg(1:4)), 1);

figure; hold on;
plot(log(res(1:4)), log(E_1_avg(1:4)), 'o-');
plot(log(res(1:4)), log(E_2_avg(1:4)), 'o-');
% plot(log(res(1:4)), log(E_3_avg(1:4)), 'o-');
% plot(log(res(1:3)), log(E_old_avg(1:3)), 'o-');
plot(log(res(1:4)), (polyval(p1, log(res(1:4)))), 'k--');
plot(log(res(1:4)), (polyval(p2, log(res(1:4)))), 'r--');
text(log(res(2)) * 0.9, log(E_1_avg(2)), ['Slope: ' num2str(p1(1)) '         '], 'HorizontalAlignment', 'right');
text(log(res(2)) * 1.1, log(E_2_avg(2)), ['         Slope: ' num2str(p2(1))], 'HorizontalAlignment', 'left');
ylabel('$L_2$ error');
xlabel('$N_{particles}$');
legend('New BC', 'New BC 2x refined wall', 'Location', 'best');

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
        PolarPosition = PolarPosition(TC.Type == TC.TypeFluid, :);
        PolarVelocity = PolarVelocity(TC.Type == TC.TypeFluid, :);
        PolarVelocityT = PolarVelocityT(TC.Type == TC.TypeFluid, :);

        [ur, uth] = TC_Analytical(PolarPosition(:, 1), [Rin Rout win wout]);

        errors(k, 1) = TC_Error(PolarPosition, PolarVelocity, ur, uth); % error in velocity
        errors(k, 2) = TC_Error(PolarPosition, PolarVelocityT, ur, uth); % error in transport velocity

        if (plotsetting == 1)
            % fine grid for plotting analytical solution plot
            rfine = linspace(Rin, Rout, 1000);
            [urfine, uthfine] = TC_Analytical(rfine, [Rin Rout win wout]);

            figure(fig1);
            sgtitle('Velocity');
            subplot(4, nres, k); hold on;
            plot(PolarPosition(:, 1), PolarVelocity(:, 2), '.'); ylabel('$u_\theta$');
            plot(rfine, uthfine);
            title([TC.PlotName "\\ $L_2$ error: " + num2str(errors(k, 1)) "$L_2$ error (transport): " + num2str(errors(k, 2))]);
            axis([1.5 2 -0.15 0.15]);

            subplot(4, nres, k + nres); hold on;
            plot(PolarPosition(:, 1), PolarVelocity(:, 1), '.'); ylabel('$u_r$');
            plot(rfine, urfine);
            axis([1.5 2 -0.1 0.1]);

            subplot(4, nres, k + 2 * nres); hold on;
            plot(PolarPosition(:, 1), abs(PolarVelocity(:, 2) - uth), '.'); ylabel('$E_\theta$');
            xlabel('$r$');
            axis([1.5 2 0 0.1]);
            title([TC.PlotName]);

            subplot(4, nres, k + 3 * nres); hold on;
            plot(PolarPosition(:, 1), abs(PolarVelocity(:, 1) - ur), '.'); ylabel('$E_r$');
            xlabel('$r$');
            axis([1.5 2 0 0.1]);

            figure(fig2);
            sgtitle('Transport velocity');
            subplot(4, nres, k); hold on;
            plot(PolarPosition(:, 1), PolarVelocityT(:, 2), '.'); ylabel('$u_\theta$');
            plot(rfine, uthfine);
            title([TC.PlotName "\\ $L_2$ error:" + num2str(errors(k, 1)) "$L_2$ error (transport): " + num2str(errors(k, 2))]);
            axis([1.5 2 -0.15 0.15]);

            subplot(4, nres, k + nres); hold on;
            plot(PolarPosition(:, 1), PolarVelocityT(:, 1), '.'); ylabel('$u_r$');
            plot(rfine, urfine);
            axis([1.5 2 -0.1 0.1]);

            subplot(4, nres, k + 2 * nres); hold on;
            plot(PolarPosition(:, 1), abs(PolarVelocityT(:, 2) - uth), '.'); ylabel('$E_\theta$');
            xlabel('$r$');
            axis([1.5 2 0 0.1]);
            title([TC.PlotName]);

            subplot(4, nres, k + 3 * nres); hold on;
            plot(PolarPosition(:, 1), abs(PolarVelocityT(:, 1) - ur), '.'); ylabel('$E_r$');
            xlabel('$r$');
            axis([1.5 2 0 0.1]);

        end

    end

    % pause
    % errors
    errors = errors(:, 2);
end
