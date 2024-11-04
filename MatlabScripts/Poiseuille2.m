clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 10);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 600 600]);

Hconst = 1;
dim = 2;
rho0 = 1;

steps = [99:99];
res = [30 60 90 120 180 240];

E_New = zeros(steps(end) - steps(1) + 1, length(res));
E_NewRef = zeros(steps(end) - steps(1) + 1, length(res));
E_Old = zeros(steps(end) - steps(1) + 1, length(res));

h = waitbar(0, 'Please wait...');
nn = length(steps);

g = 0.1;
nu = 0.1;
L = 1;
params = [g nu L];

umax_new = zeros(steps(end) - steps(1) + 1, length(res));
umax_newref = zeros(steps(end) - steps(1) + 1, length(res));
umax_old = zeros(steps(end) - steps(1) + 1, length(res));

%%%%%%%%%%%% FOR LOADING NEW DATA %%%%%%%%%%%%%
for k = 1:length(steps)
    nfile = steps(k);

    dataNew1 = ParticleData('../CSV_Data/Poiseuille_new_summation_12_30_1rf_1prc/file', nfile, ['30 new'], rho0, dim, Hconst);
    dataNew2 = ParticleData('../CSV_Data/Poiseuille_new_summation_24_60_1rf_1prc/file', nfile, ['60 new'], rho0, dim, Hconst);
    dataNew3 = ParticleData('../CSV_Data/Poiseuille_new_summation_36_90_1rf_1prc/file', nfile, ['90 new'], rho0, dim, Hconst);
    dataNew4 = ParticleData('../CSV_Data/Poiseuille_new_summation_48_120_1rf_1prc/file', nfile, ['120 new'], rho0, dim, Hconst);
    dataNew5 = ParticleData('../CSV_Data/Poiseuille_new_summation_72_180_1rf_1prc/file', nfile, ['180 new'], rho0, dim, Hconst);
    dataNew6 = ParticleData('../CSV_Data/Poiseuille_new_summation_96_240_1rf_1prc/file', nfile, ['240 new'], rho0, dim, Hconst);

    dataNewRef1 = ParticleData('../CSV_Data/Poiseuille_new_summation_12_30_2rf_1prc/file', nfile, ['30 new refined'], rho0, dim, Hconst);
    dataNewRef2 = ParticleData('../CSV_Data/Poiseuille_new_summation_24_60_2rf_1prc/file', nfile, ['60 new refined'], rho0, dim, Hconst);
    dataNewRef3 = ParticleData('../CSV_Data/Poiseuille_new_summation_36_90_2rf_1prc/file', nfile, ['90 new refined'], rho0, dim, Hconst);
    dataNewRef4 = ParticleData('../CSV_Data/Poiseuille_new_summation_48_120_2rf_1prc/file', nfile, ['120 new refined'], rho0, dim, Hconst);
    dataNewRef5 = ParticleData('../CSV_Data/Poiseuille_new_summation_72_180_2rf_1prc/file', nfile, ['180 new refined'], rho0, dim, Hconst);
    dataNewRef6 = ParticleData('../CSV_Data/Poiseuille_new_summation_96_240_2rf_1prc/file', nfile, ['240 new refined'], rho0, dim, Hconst);

    dataOld1 = ParticleData('../CSV_Data/Poiseuille_old_summation_12_30_1rf_1prc/file', nfile, ['30 old'], rho0, dim, Hconst);
    dataOld2 = ParticleData('../CSV_Data/Poiseuille_old_summation_24_60_1rf_1prc/file', nfile, ['60 old'], rho0, dim, Hconst);
    dataOld3 = ParticleData('../CSV_Data/Poiseuille_old_summation_36_90_1rf_1prc/file', nfile, ['90 old'], rho0, dim, Hconst);
    dataOld4 = ParticleData('../CSV_Data/Poiseuille_old_summation_48_120_1rf_1prc/file', nfile, ['120 old'], rho0, dim, Hconst);
    dataOld5 = ParticleData('../CSV_Data/Poiseuille_old_summation_72_180_1rf_1prc/file', nfile, ['180 old'], rho0, dim, Hconst);
    dataOld6 = ParticleData('../CSV_Data/Poiseuille_old_summation_96_240_1rf_1prc/file', nfile, ['240 old'], rho0, dim, Hconst);

    setNew = {dataNew1 dataNew2 dataNew3 dataNew4 dataNew5 dataNew6};
    setNewRef = {dataNewRef1 dataNewRef2 dataNewRef3 dataNewRef4 dataNewRef5 dataNewRef6};
    setOld = {dataOld1 dataOld2 dataOld3 dataOld4 dataOld5 dataOld6};

    umax_new(k, :) = [max(dataNew1.Velocity(:, 1)) max(dataNew2.Velocity(:, 1)) max(dataNew3.Velocity(:, 1)) max(dataNew4.Velocity(:, 1)) max(dataNew5.Velocity(:, 1)) max(dataNew6.Velocity(:, 1))];
    umax_newref(k, :) = [max(dataNewRef1.Velocity(:, 1)) max(dataNewRef2.Velocity(:, 1)) max(dataNewRef3.Velocity(:, 1)) max(dataNewRef4.Velocity(:, 1)) max(dataNewRef5.Velocity(:, 1)) max(dataNewRef6.Velocity(:, 1))];
    umax_old(k, :) = [max(dataOld1.Velocity(:, 1)) max(dataOld2.Velocity(:, 1)) max(dataOld3.Velocity(:, 1)) max(dataOld4.Velocity(:, 1)) max(dataOld5.Velocity(:, 1)) max(dataOld6.Velocity(:, 1))];

    plot_setting = 0;

    if (k == length(steps))
        plot_setting = 1;
    end

    E_New(k, :) = ComputeErrorAndPlotPoiseuille(setNew, plot_setting, params);
    E_NewRef(k, :) = ComputeErrorAndPlotPoiseuille(setNewRef, plot_setting, params);
    E_Old(k, :) = ComputeErrorAndPlotPoiseuille(setOld, plot_setting, params);

    waitbar(k / nn, h, sprintf('Progress: %d percent ..', round((k / nn) * 100)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('Poiseuille2.mat', 'E_New', 'E_Old', 'umax_new', 'umax_old', 'res', 'steps');
% load('Poiseuille2.mat');
close(h);

figure; hold on;
plot(steps, umax_new(:, 1), 'r-', 'DisplayName', '30 new');
plot(steps, umax_new(:, 2), 'g-', 'DisplayName', '60 new');
plot(steps, umax_new(:, 3), 'b-', 'DisplayName', '90 new');
plot(steps, umax_new(:, 4), 'k-', 'DisplayName', '120 new');
plot(steps, umax_new(:, 5), 'm-', 'DisplayName', '180 new');
plot(steps, umax_new(:, 6), 'y-', 'DisplayName', '240 new');
% plot(steps, umax_old(:, 1), 'r:', 'DisplayName', '30 old');
% plot(steps, umax_old(:, 2), 'g:', 'DisplayName', '60 old');
% plot(steps, umax_old(:, 3), 'b:', 'DisplayName', '90 old');
% plot(steps, umax_old(:, 4), 'k:', 'DisplayName', '120 old');
% plot(steps, umax_old(:, 5), 'm:', 'DisplayName', '180 old');
% plot(steps, umax_old(:, 6), 'y:', 'DisplayName', '240 old');
yline(0.125, 'k--', 'DisplayName', 'Analytical');
legend('Location', 'best');

figure; hold on;
plot(steps, log(E_New(:, 1)), 'r-', 'DisplayName', '30 new');
plot(steps, log(E_New(:, 2)), 'g-', 'DisplayName', '60 new');
plot(steps, log(E_New(:, 3)), 'b-', 'DisplayName', '90 new');
plot(steps, log(E_New(:, 4)), 'k-', 'DisplayName', '120 new');
plot(steps, log(E_New(:, 5)), 'm-', 'DisplayName', '180 new');
plot(steps, log(E_New(:, 6)), 'y-', 'DisplayName', '240 new');
% plot(steps, log(E_Old(:, 1)), 'r:', 'DisplayName', '30 old');
% plot(steps, log(E_Old(:, 2)), 'g:', 'DisplayName', '60 old');
% plot(steps, log(E_Old(:, 3)), 'b:', 'DisplayName', '90 old');
% plot(steps, log(E_Old(:, 4)), 'k:', 'DisplayName', '120 old');
% plot(steps, log(E_Old(:, 5)), 'm:', 'DisplayName', '180 old');
% plot(steps, log(E_Old(:, 6)), 'y:', 'DisplayName', '240 old');
legend('Location', 'best');

% central_step = 100; avg_window = 2;
% E_New_avg = mean(E_New(central_step - avg_window:central_step + avg_window, :), 1);
% E_Old_avg = mean(E_Old(central_step - avg_window:central_step + avg_window, :), 1);

E_New_avg = E_New(end, :);
E_NewRef_avg = E_NewRef(end, :);
E_Old_avg = E_Old(end, :);

pnew = polyfit(log(res), log(E_New_avg), 1);
pnewref = polyfit(log(res), log(E_NewRef_avg), 1);
pold = polyfit(log(res), log(E_Old_avg), 1);

figure; hold on;
plot(log(res), log(E_New_avg), 'o-');
plot(log(res), log(E_NewRef_avg), 'o-');
plot(log(res), log(E_Old_avg), 'o-');
plot(log(res), pnew(1) * log(res) + pnew(2), 'k--');
plot(log(res), pnewref(1) * log(res) + pnewref(2), 'k--');
plot(log(res), pold(1) * log(res) + pold(2), 'k--');

% Define the range for the reference lines based on log(res)
x_range = [min(log(res)), max(log(res))];

% Set offsets to adjust the position of the reference lines
offset_1 = 0.5; % Adjust this as needed
offset_2 = -0.5; % Adjust this as needed

% Reference line with slope -1, shifted up
y1 = log(E_New_avg(1)) + (-1) * (x_range - log(res(1))) + offset_1;
plot(x_range, y1, 'r');

% Reference line with slope -2, shifted down
y2 = log(E_New_avg(1)) + (-2) * (x_range - log(res(1))) + offset_2;
plot(x_range, y2, 'r:');

legend('E New', 'E Old', ['Enew fit slope =', num2str(pnew(1))], ['Eold fit slope = ', num2str(pold(1))], 'Slope -1', 'Slope -2');
xlabel('$\log(N_y)$');
ylabel('log(L_1 error)');
title('Error Convergence with Reference Slopes');
hold off;

function [ux, uy] = Poiseuille_Analytical(y, params)
    g = params(1);
    nu = params(2);
    L = params(3);
    ux = (g / (2 * nu)) .* y .* (L - y);
    uy = zeros(size(y));
end

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

function [L2, L2spatial] = L2_error(VelocityNumeric, VelocityAnalytic)
    L2spatial = (VelocityNumeric(:, 1) - VelocityAnalytic(:, 1)) .^ 2 + (VelocityNumeric(:, 2) - VelocityAnalytic(:, 2)) .^ 2;
    N = length(VelocityNumeric(:, 1));
    L2 = sqrt((1 / N) * sum(L2spatial));
end

function [L2, L2spatial] = L2x_error(VelocityNumeric, VelocityAnalytic)
    L2spatial = (VelocityNumeric(:, 1) - VelocityAnalytic(:, 1)) .^ 2;
    N = length(VelocityNumeric(:, 1));
    L2 = sqrt((1 / N) * sum(L2spatial));
end

function [L2, L2spatial] = L2y_error(VelocityNumeric, VelocityAnalytic)
    L2spatial = (VelocityNumeric(:, 2) - VelocityAnalytic(:, 2)) .^ 2;
    N = length(VelocityNumeric(:, 2));
    L2 = sqrt((1 / N) * sum(L2spatial));
end

function [L1, L1spatial] = L1_error(VelocityNumeric, VelocityAnalytic)
    L1spatial = abs(VelocityNumeric(:, 1) - VelocityAnalytic(:, 1)) + abs(VelocityNumeric(:, 2) - VelocityAnalytic(:, 2));
    N = length(VelocityNumeric(:, 1));
    L1 = (1 / N) * sum(L1spatial);
end

function [L1, L1spatial] = L1x_error(VelocityNumeric, VelocityAnalytic)
    L1spatial = abs(VelocityNumeric(:, 1) - VelocityAnalytic(:, 1));
    N = length(VelocityNumeric(:, 1));
    L1 = (1 / N) * sum(L1spatial);
end

function [L1, L1spatial] = L1y_error(VelocityNumeric, VelocityAnalytic)
    L1spatial = abs(VelocityNumeric(:, 2) - VelocityAnalytic(:, 2));
    N = length(VelocityNumeric(:, 2));
    L1 = (1 / N) * sum(L1spatial);
end

function [Erel, ErelSpatial] = Relative_error(VelocityNumeric, VelocityAnalytic)
    ErelSpatial = sqrt((VelocityNumeric(:, 1) - VelocityAnalytic(:, 1)) .^ 2 + (VelocityNumeric(:, 2) - VelocityAnalytic(:, 2)) .^ 2) ./ sqrt(VelocityAnalytic(:, 1) .^ 2 + VelocityAnalytic(:, 2) .^ 2);
    N = length(VelocityNumeric(:, 1));
    Erel = (1 / N) * sum(ErelSpatial);
end

function [Erel, ErelSpatial] = Relative_errorX(VelocityNumeric, VelocityAnalytic)
    ErelSpatial = sqrt((VelocityNumeric(:, 1) - VelocityAnalytic(:, 1)) .^ 2) ./ sqrt(VelocityAnalytic(:, 1) .^ 2);
    N = length(VelocityNumeric(:, 1));
    Erel = (1 / N) * sum(ErelSpatial);
end

% function [MSE, MSEx] = MSE_error(VelocityNumeric, VelocityAnalytic)
%     MSEx = (VelocityNumeric(:, 1) - VelocityAnalytic(:, 1)) .^ 2 + (VelocityNumeric(:, 2) - VelocityAnalytic(:, 2)) .^ 2;
%     N = length(VelocityNumeric(:, 1));
%     MSE = (1 / N) * sum(MSEx);
% end

% function [L2, spat] = P_Error(Velocity, ux, uy)
%     NormSPH = sqrt(Velocity(:, 1) .^ 2 + Velocity(:, 2) .^ 2);
%     NormAnalytical = sqrt(ux .^ 2 + uy .^ 2);
%     spat = (NormAnalytical - NormSPH) .^ 2;

%     L2 = sqrt((1 / length(Velocity(:, 1))) * sum (spat));
% end

function [errors] = ComputeErrorAndPlotPoiseuille(datasets, plot_setting, params)

    nres = length(datasets); % number of resolutions in the dataset
    errors = zeros(nres, 2);

    g = params(1); nu = params(2); L = params(3);

    for k = 1:nres
        data = datasets{k};

        pos = data.Position;
        vel = data.Velocity;
        velT = data.VelocityTransport;

        % Remove boundary particles
        pos = pos(data.Type == data.TypeFluid, :);
        vel = vel(data.Type == data.TypeFluid, :);
        velT = velT(data.Type == data.TypeFluid, :);

        % Remove z component (all zeros)
        pos = pos(:, 1:2);
        vel = vel(:, 1:2);
        velT = velT(:, 1:2);

        % sort the particles based on y
        [pos, I] = sortrows(pos, 2);
        vel = vel(I, :);
        velT = velT(I, :);

        % get H
        H = data.H;

        % restrict data to y < 3H and y > L - 3H
        % subset_idx = pos(:, 2) < 3 * H | pos(:, 2) > L - 3 * H;
        % pos = pos(subset_idx, :);
        % vel = vel(subset_idx, :);
        % velT = velT(subset_idx, :);

        [ux, uy] = Poiseuille_Analytical(pos(:, 2), params);
        velA = [ux uy];

        [ErrorMagnitude, ErrorSpatial] = L1x_error(vel, velA);
        [ErrorMagnitudeT, ErrorSpatialT] = L1x_error(velT, velA);

        errors(k, 1) = ErrorMagnitude;
        errors(k, 2) = ErrorMagnitudeT;

        if (plot_setting == 0)
            continue;
        end

        % fine grid for analytical solution plot
        yfine = linspace(0, 1, 1000);
        [uxfine, uyfine] = Poiseuille_Analytical(yfine, params);

        figure; sgtitle(data.PlotName);
        subplot(3, 2, 1); hold on; % plot u(y)
        plot(yfine, uxfine, 'k-', 'DisplayName', 'Analytical');
        plot(pos(:, 2), vel(:, 1), '-ob', 'DisplayName', 'SPH');
        xlabel('$y$'); ylabel('$u_x$');
        legend('Location', 'best');
        title('Momentum velocity');

        subplot(3, 2, 3); hold on; % plot v(y)
        plot(yfine, uyfine, 'k-', 'DisplayName', 'Analytical');
        plot(pos(:, 2), vel(:, 2), '-ob', 'DisplayName', 'SPH');
        xlabel('$y$'); ylabel('$u_y$');
        legend('Location', 'best');

        subplot(3, 2, 5); hold on; % plot error(y)
        plot(pos(:, 2), ErrorSpatial, '-ob', 'DisplayName', 'L_1 error');
        xlabel('$y$'); ylabel('L1 x');

        %% transport velocity
        subplot(3, 2, 2); hold on; % plot u(y)
        plot(yfine, uxfine, 'k-', 'DisplayName', 'Analytical');
        plot(pos(:, 2), velT(:, 1), '-ob', 'DisplayName', 'SPH');
        xlabel('$y$'); ylabel('$u_x$');
        legend('Location', 'best');
        title('Transport velocity');

        subplot(3, 2, 4); hold on; % plot v(y)
        plot(yfine, uyfine, 'k-', 'DisplayName', 'Analytical');
        plot(pos(:, 2), velT(:, 2), '-ob', 'DisplayName', 'SPH');
        xlabel('$y$'); ylabel('$u_y$');
        legend('Location', 'best');

        subplot(3, 2, 6); hold on; % plot error(y)
        plot(pos(:, 2), ErrorSpatialT, '-ob', 'DisplayName', 'Relative error');
        xlabel('$y$'); ylabel('L1 x');

    end

    errors = errors(:, 1);
end
