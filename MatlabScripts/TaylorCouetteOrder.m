clc; close all;
global TypeBoundary TypeObstacle TypeFluid

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultLineMarkerSize', 2);
set(groot, 'defaultFigurePosition', [1440 0 600 600]);

Hconst = 1; dim = 2; rho0 = 1;

steps = [300:300];
res = [10 20 30 40 50];
hsize = 0.5 ./ res;
dirname = 'TCNEW';

E_1 = zeros(steps(end) - steps(1) + 1, length(res));
E_2 = zeros(steps(end) - steps(1) + 1, length(res));
% E_05 = zeros(steps(end) - steps(1) + 1, length(res));
% E_sum = zeros(steps(end) - steps(1) + 1, length(res));
E_old = zeros(steps(end) - steps(1) + 1, length(res));

Rin = 0.5;
Rout = 1.0;
Win = 1.0;
Wout = 1.0;
params = [Rin, Rout, Win, Wout];

nn = length(steps);
h = waitbar(0, 'Please wait...');

for k = 1:length(steps)
    nfile = steps(k);

    tc_10_1 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_10_10_1rf_1prc/file'], nfile, ['10x1'], rho0, dim, Hconst);
    tc_10_2 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_10_10_2rf_1prc/file'], nfile, ['10x2'], rho0, dim, Hconst);
    % tc_10_05 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_flatvolume_10_10_1rf_1prc/file'], nfile, ['10x0.5'], rho0, dim, Hconst);
    % tc_10_sum = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_10_10_1rf_1prc/file'], nfile, ['10 sum'], rho0, dim, Hconst);
    tc_10_old = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_old_differential_10_10_1rf_1prc/file'], nfile, ['10 old'], rho0, dim, Hconst);

    tc_20_1 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_20_20_1rf_1prc/file'], nfile, ['20x1'], rho0, dim, Hconst);
    tc_20_2 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_20_20_2rf_1prc/file'], nfile, ['20x2'], rho0, dim, Hconst);
    % tc_20_05 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_flatvolume_20_20_1rf_1prc/file'], nfile, ['20x0.5'], rho0, dim, Hconst);
    % tc_20_sum = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_20_20_1rf_1prc/file'], nfile, ['20 sum'], rho0, dim, Hconst);
    tc_20_old = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_old_differential_20_20_1rf_1prc/file'], nfile, ['20 old'], rho0, dim, Hconst);

    tc_30_1 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_30_30_1rf_1prc/file'], nfile, ['30x1'], rho0, dim, Hconst);
    tc_30_2 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_30_30_2rf_1prc/file'], nfile, ['30x2'], rho0, dim, Hconst);
    % tc_30_05 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_flatvolume_30_30_1rf_1prc/file'], nfile, ['30x0.5'], rho0, dim, Hconst);
    % tc_30_sum = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_30_30_1rf_1prc/file'], nfile, ['30 sum'], rho0, dim, Hconst);
    tc_30_old = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_old_differential_30_30_1rf_1prc/file'], nfile, ['30 old'], rho0, dim, Hconst);

    tc_40_1 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_40_40_1rf_1prc/file'], nfile, ['40x1'], rho0, dim, Hconst);
    tc_40_2 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_40_40_2rf_1prc/file'], nfile, ['40x2'], rho0, dim, Hconst);
    % tc_40_05 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_flatvolume_40_40_1rf_1prc/file'], nfile, ['40x0.5'], rho0, dim, Hconst);
    % tc_40_sum = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_40_40_1rf_1prc/file'], nfile, ['40 sum'], rho0, dim, Hconst);
    tc_40_old = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_old_differential_40_40_1rf_1prc/file'], nfile, ['40 old'], rho0, dim, Hconst);

    tc_50_1 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_50_50_1rf_1prc/file'], nfile, ['50x1'], rho0, dim, Hconst);
    tc_50_2 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_50_50_2rf_1prc/file'], nfile, ['50x2'], rho0, dim, Hconst);
    % tc_50_05 = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_differential_flatvolume_50_50_1rf_1prc/file'], nfile, ['50x0.5'], rho0, dim, Hconst);
    % tc_50_sum = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_new_summation_50_50_1rf_1prc/file'], nfile, ['50 sum'], rho0, dim, Hconst);
    tc_50_old = ParticleData(['../CSV_Data/' dirname '/TaylorCouette_old_differential_50_50_1rf_1prc/file'], nfile, ['50 old'], rho0, dim, Hconst);

    set_1 = {tc_10_1, tc_20_1, tc_30_1 tc_40_1, tc_50_1};
    set_2 = {tc_10_2 tc_20_2, tc_30_2 tc_40_2, tc_50_2};
    % set_05 = {tc_10_05, tc_20_05, tc_30_05 tc_40_05, tc_50_05};
    % set_sum = {tc_10_sum, tc_20_sum, tc_30_sum tc_40_sum, tc_50_sum};
    set_old = {tc_10_old, tc_20_old, tc_30_old tc_40_old, tc_50_old};

    plotsetting = 0;
    restrictionsetting = 0;

    % if k == length(steps)
    %     plotsetting = 1;
    % end

    % Error of the 10, 20, 30 resolutions at time step nfile

    E_1(k, :) = ComputeErrorAndPlotTC(set_1, params, plotsetting, restrictionsetting);
    E_2(k, :) = ComputeErrorAndPlotTC(set_2, params, plotsetting, restrictionsetting);
    % E_05(k, :) = ComputeErrorAndPlotTC(set_05, params, plotsetting, restrictionsetting);
    % E_sum(k, :) = ComputeErrorAndPlotTC(set_sum, params, plotsetting, restrictionsetting);
    E_old(k, :) = ComputeErrorAndPlotTC(set_old, params, plotsetting, restrictionsetting);

    disp(['Step: ' num2str(nfile)]);
    waitbar(k / nn, h, sprintf('Progress: %d percent ..', round((k / nn) * 100)));

end

close(h);

% save('TCOrder.mat', 'E_1', 'E_2', 'E_old', 'res', 'steps');
% load('TCOrder.mat');
% figure; hold on;
% plot(steps, log10(E_1(:, 1)), 'o-');
% plot(steps, log10(E_1(:, 2)), 'o-');
% plot(steps, log10(E_1(:, 3)), 'o-');
% plot(steps, log10(E_1(:, 4)), 'o-');
% plot(steps, log10(E_1(:, 5)), 'o-');
% title('x1'); legend('10', '20', '30', '40', '50');

% figure; hold on;
% plot(steps, log10(E_2(:, 1)), 'o-');
% plot(steps, log10(E_2(:, 2)), 'o-');
% plot(steps, log10(E_2(:, 3)), 'o-');
% plot(steps, log10(E_2(:, 4)), 'o-');
% plot(steps, log10(E_2(:, 5)), 'o-');
% title('x2'); legend('10', '20', '30', '40', '50');

% figure; hold on;
% plot(steps, log10(E_old(:, 1)), 'o-');
% plot(steps, log10(E_old(:, 2)), 'o-');
% plot(steps, log10(E_old(:, 3)), 'o-');
% plot(steps, log10(E_old(:, 4)), 'o-');
% plot(steps, log10(E_old(:, 5)), 'o-');
% title('old'); legend('10', '20', '30', '40', '50');

% figure; hold on;
% plot(steps, log10(E_sum(:, 1)), 'o-');
% plot(steps, log10(E_sum(:, 2)), 'o-');
% plot(steps, log10(E_sum(:, 3)), 'o-');
% plot(steps, log10(E_sum(:, 4)), 'o-');
% plot(steps, log10(E_sum(:, 5)), 'o-');
% title('sum'); legend('10', '20', '30', '40', '50');

E_1_avg = E_1(end, :);
E_2_avg = E_2(end, :);
E_05_avg = E_05(end, :);
% E_sum_avg = E_sum(end, :);
E_old_avg = E_old(end, :);

p1 = polyfit(log10(hsize), log10(E_1_avg), 1);
p2 = polyfit(log10(hsize), log10(E_2_avg), 1);
% p05 = polyfit(log10(hsize), log10(E_05_avg), 1);
% psum = polyfit(log10(hsize), log10(E_sum_avg), 1);
pold = polyfit(log10(hsize), log10(E_old_avg), 1);

% Define the range for the reference lines based on log10(hsize)
x_range = [min(log10(hsize)), max(log10(hsize))];
offset_1 = 0.5;
offset_2 = -0.5;
% Reference line with slope 1, shifted up
y1 = log10(E_old_avg(1)) + (1) * (x_range - log10(hsize(1))) + offset_1;
% Reference line with slope 2, shifted down
y2 = log10(E_1_avg(1)) + (2) * (x_range - log10(hsize(1))) + offset_2;

figure; hold on;
plot(log10(hsize), log10(E_1_avg), '-s', 'MarkerSize', 10);
plot(log10(hsize), log10(E_2_avg), '-s', 'MarkerSize', 10);
% plot(log10(hsize), log10(E_05_avg), '-o');
% plot(log10(hsize), log10(E_sum_avg), '-o');
plot(log10(hsize), log10(E_old_avg), '-s', 'MarkerSize', 10);
% plot(log10(res), p1(1) * log10(res) + p1(2), 'r--');
% plot(log10(res), p2(1) * log10(res) + p2(2), 'b--');
% plot(log10(res), pold(1) * log10(res) + pold(2), 'k--');
plot(x_range, y1, 'r');
plot(x_range, y2, 'k');
ylabel('$\log_{10}(L_1)$');
xlabel('$\log_{10}(h)$');
legend('New BC', 'New BC 2x', 'Old BC', 'slope 1', 'slope 2', 'Location', 'bestoutside');
% legend('New BC', 'Old BC', 'slope -1', 'slope -2');

% put text in lower right corner of the plot
pos = get(gca, 'Position');
txt = {['slope new x1 = ' num2str(p1(1))]; ['slope new x2 = ' num2str(p2(1))]; [' slope old = ' num2str(pold(1))]};
text(pos(1) + 0.6, pos(2) + 0.1, txt, 'Units', 'normalized');

function [ur, uth] = TC_Analytical(r, params)
    Rin = params(1);
    Rout = params(2);
    Win = params(3);
    Wout = params(4);
    a =- ((Rout ^ 2 * Rin ^ 2) / (Rout ^ 2 - Rin ^ 2)) * (Wout - Win);
    b = (Wout * Rout ^ 2 - Win * Rin ^ 2) / (Rout ^ 2 - Rin ^ 2);
    ur = zeros(length(r), 1);
    uth = a ./ r + b * r;

end

%% L2 error
function [L2, L2spatial] = L2_error(VelocityNumeric, VelocityAnalytic)
    L2spatial = (VelocityNumeric(:, 1) - VelocityAnalytic(:, 1)) .^ 2 + (VelocityNumeric(:, 2) - VelocityAnalytic(:, 2)) .^ 2;
    N = length(VelocityNumeric(:, 1));
    L2 = sqrt((1 / N) * sum(L2spatial));
end

function [L2, L2spatial] = L2r_error(VelocityNumeric, VelocityAnalytic)
    L2spatial = (VelocityNumeric(:, 1) - VelocityAnalytic(:, 1)) .^ 2;
    N = length(VelocityNumeric(:, 1));
    L2 = sqrt((1 / N) * sum(L2spatial));
end

function [L2, L2spatial] = L2th_error(VelocityNumeric, VelocityAnalytic)
    L2spatial = (VelocityNumeric(:, 2) - VelocityAnalytic(:, 2)) .^ 2;
    N = length(VelocityNumeric(:, 2));
    L2 = sqrt((1 / N) * sum(L2spatial));
end

%%

%% L1 error
function [L1, L1spatial] = L1_error(VelocityNumeric, VelocityAnalytic)
    L1spatial = abs(VelocityNumeric(:, 1) - VelocityAnalytic(:, 1)) + abs(VelocityNumeric(:, 2) - VelocityAnalytic(:, 2));
    N = length(VelocityNumeric(:, 1));
    L1 = (1 / N) * sum(L1spatial);
end

function [L1, L1spatial] = L1r_error(VelocityNumeric, VelocityAnalytic)
    L1spatial = abs(VelocityNumeric(:, 1) - VelocityAnalytic(:, 1));
    N = length(VelocityNumeric(:, 1));
    L1 = (1 / N) * sum(L1spatial);
end

function [L1, L1spatial] = L1th_error(VelocityNumeric, VelocityAnalytic)
    L1spatial = abs(VelocityNumeric(:, 2) - VelocityAnalytic(:, 2));
    N = length(VelocityNumeric(:, 2));
    L1 = (1 / N) * sum(L1spatial);
end

%%

%% Relative error
function [Erel, ErelSpatial] = Relative_error(VelocityNumeric, VelocityAnalytic)
    ErelSpatial = sqrt((VelocityNumeric(:, 1) - VelocityAnalytic(:, 1)) .^ 2 + (VelocityNumeric(:, 2) - VelocityAnalytic(:, 2)) .^ 2) ./ sqrt(VelocityAnalytic(:, 1) .^ 2 + VelocityAnalytic(:, 2) .^ 2);
    N = length(VelocityNumeric(:, 1));
    Erel = (1 / N) * sum(ErelSpatial);
end

function [Erel, ErelSpatial] = Relative_errorR(VelocityNumeric, VelocityAnalytic)
    ErelSpatial = sqrt((VelocityNumeric(:, 1) - VelocityAnalytic(:, 1)) .^ 2) ./ sqrt(VelocityAnalytic(:, 1) .^ 2);
    N = length(VelocityNumeric(:, 1));
    Erel = (1 / N) * sum(ErelSpatial);
end

%%

% function [L2] = TC_Error(PolarPosition, PolarVelocity, ur, uth)
%     NormSPH = sqrt(PolarVelocity(:, 1) .^ 2 + PolarVelocity(:, 2) .^ 2);
%     NormAnalytical = sqrt(ur .^ 2 + uth .^ 2);

%     L2 = sqrt((1 / length(PolarPosition(:, 1))) * sum ((NormAnalytical - NormSPH) .^ 2));
% end

% function [L2] = R_Error(PolarPosition, PolarVelocity, ur, uth)
%     NormSPH = sqrt(PolarVelocity(:, 1) .^ 2);
%     NormAnalytical = sqrt(ur .^ 2);

%     L2 = sqrt((1 / length(PolarPosition(:, 1))) * sum ((NormAnalytical - NormSPH) .^ 2));
% end

% function [L2] = TH_Error(PolarPosition, PolarVelocity, ur, uth)
%     NormSPH = sqrt(PolarVelocity(:, 2) .^ 2);
%     NormAnalytical = sqrt(uth .^ 2);

%     L2 = sqrt((1 / length(PolarPosition(:, 1))) * sum ((NormAnalytical - NormSPH) .^ 2));
% end

function [errors] = ComputeErrorAndPlotTC(datasets, params, plotsetting, restriction_setting)

    nres = length(datasets);
    errors = zeros(nres, 2);

    Rin = params(1); Rout = params(2); Win = params(3); Wout = params(4);

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

        % sort based on r
        [PolarPosition, I] = sortrows(PolarPosition, 1);
        PolarVelocity = PolarVelocity(I, :);
        PolarVelocityT = PolarVelocityT(I, :);

        %  restrict data to r < Rin+ 3H and R > Rout - 3H
        if (restriction_setting == 1)
            H = TC.H;
            subset_idx = PolarPosition(:, 1) < Rin + 3 * H | PolarPosition(:, 1) > Rout - 3 * H;
            PolarPosition = PolarPosition(subset_idx, :);
            PolarVelocity = PolarVelocity(subset_idx, :);
            PolarVelocityT = PolarVelocityT(subset_idx, :);
        end

        [ur, uth] = TC_Analytical(PolarPosition(:, 1), params);
        PolarVelA = [ur, uth];

        [L1, L1spatial] = L1_error(PolarVelocity, PolarVelA); % error in velocity
        [L1T, L1spatialT] = L1_error(PolarVelocityT, PolarVelA); % error in velocity transport

        errors(k, 1) = L1;
        errors(k, 2) = L1T;

        if (plotsetting == 0)
            continue;
        end

        % fine grid for plotting analytical solution plot
        rfine = linspace(Rin, Rout, 1000);
        [urfine, uthfine] = TC_Analytical(rfine, params);

        figure; sgtitle(TC.PlotName);
        subplot(3, 2, 1); hold on;
        plot(rfine, uthfine, 'k-', 'DisplayName', 'Analytical');
        plot(PolarPosition(:, 1), PolarVelocity(:, 2), 'ob', 'DisplayName', 'SPH');
        xlabel('$r$'); ylabel('$u_\theta$');
        legend('location', 'best');
        title('Momentum Velocity');

        subplot(3, 2, 3); hold on;
        plot(rfine, urfine, 'k-', 'DisplayName', 'Analytical');
        plot(PolarPosition(:, 1), PolarVelocity(:, 1), '.b', 'DisplayName', 'SPH');
        xlabel('$r$'); ylabel('$u_r$');
        legend('location', 'best');

        subplot(3, 2, 5); hold on;
        plot(PolarPosition(:, 1), L1spatial, '-.b');
        xlabel('$r$'); ylabel('$L_1$');

        % transport velocity
        subplot(3, 2, 2); hold on;
        plot(rfine, uthfine, 'k-', 'DisplayName', 'Analytical');
        plot(PolarPosition(:, 1), PolarVelocityT(:, 2), 'ob', 'DisplayName', 'SPH');
        xlabel('$r$'); ylabel('$u_\theta$');
        legend('location', 'best');
        title('Transport Velocity');

        subplot(3, 2, 4); hold on;
        plot(rfine, urfine, 'k-', 'DisplayName', 'Analytical');
        plot(PolarPosition(:, 1), PolarVelocityT(:, 1), '.b', 'DisplayName', 'SPH');
        xlabel('$r$'); ylabel('$u_r$');
        legend('location', 'best');

        subplot(3, 2, 6); hold on;
        plot(PolarPosition(:, 1), L1spatialT, '-.b');
        xlabel('$r$'); ylabel('$L_1$');

    end

    errors = errors(:, 1);
end
