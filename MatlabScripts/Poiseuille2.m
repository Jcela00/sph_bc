clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 10);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 600 600]);
dp = 10000;
Nfluid = [20 40 1];
Nboundary = [0, 0, 0];
Hconst = 1;
dim = 2;
rho0 = 1;

steps = 1:250;

% E_1 = zeros(length(steps), 2);
% E_2 = zeros(length(steps), 2);
% E_3 = zeros(length(steps), 2);
% E_4 = zeros(length(steps), 2);
% E_old = zeros(length(steps), 2);

% for k = 1:length(steps)
%     k
nfile = 250;
P_20_40_1 = ParticleData('../CSV_Data/Poiseuille_NEW_BC_Summation_20_40_4prc_1rf/file', nfile, ['20x40x1'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
P_20_40_old = ParticleData('../CSV_Data/Poiseuille_OLD_BC_Summation_20_40_4prc/file', nfile, ['20x40 old'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

P_40_80_1 = ParticleData('../CSV_Data/Poiseuille_NEW_BC_Summation_40_80_4prc_1rf/file', nfile, ['40x80x1'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
P_40_80_old = ParticleData('../CSV_Data/Poiseuille_OLD_BC_Summation_40_80_4prc/file', nfile, ['40x80 old'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

% P_60_120_1 = ParticleData('../CSV_Data/Poiseuille_NEW_BC_Summation_60_120_4prc_1rf/file', nfile, ['60x120x1'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
% P_60_120_old = ParticleData('../CSV_Data/Poiseuille_OLD_BC_Summation_60_120_4prc/file', nfile, ['60x120 old'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

set20 = {P_20_40_1, P_20_40_old};
set40 = {P_40_80_1, P_40_80_old};
% set60 = {P_60_120_1, P_60_120_old};

E_20 = ComputeErrorAndPlotPoiseuille(set20, 0);
E_40 = ComputeErrorAndPlotPoiseuille(set40, 0);
% E_60 = ComputeErrorAndPlotPoiseuille(set60, 0);

% E_1(k, :) = [E_20(1) E_40(1)];
% E_old(k, :) = [E_20(2) E_40(2)];
E_1 = [E_20(1) E_40(1)];
E_old = [E_20(2) E_40(2)];

% % do linear fit of E_1 and E_old
res = [20 40];
p1 = polyfit(log(res), log(E_1), 1);
p2 = polyfit(log(res), log(E_old), 1);

% figure;
% loglog(res, E_1, 'o-'); hold on;
% loglog(res, E_2, 'o-');
% loglog(res, E_3, 'o-');
% loglog(res, E_4, 'o-');
% loglog(res, E_old, 'o-');
% ylabel('$L_2$ error');
% xlabel('$N_{particles}$');
% legend('New BC', 'New BC 2x refined wall', 'New BC 3x refined wall', 'New BC 4x refined wall', 'Old BC', 'Location', 'best');

figure;
loglog(res, E_1, 'o-'); hold on;
loglog(res, exp(polyval(p1, log(res))), 'k--');
loglog(res, E_old, 'o-');
loglog(res, exp(polyval(p2, log(res))), 'r--');
text(res(2) * 0.9, E_1(2), ['Slope: ' num2str(p1(1)) '         '], 'HorizontalAlignment', 'right');
text(res(2) * 1.1, E_old(2), ['         Slope: ' num2str(p2(1))], 'HorizontalAlignment', 'left');
ylabel('$L_2$ error');
xlabel('$N_{particles}$');
legend('New BC', 'New BC linear fit', 'Old BC', 'Old BC linear fit', 'Location', 'best');

% E_20 = E_20(1:end - 1);
% E_40 = E_40(1:end - 1);
% E_60 = E_60(1:end - 1);

% figure;
% ref = [1 2 3 4];
% loglog(ref, E_20, 'o-');
% % loglog(ref, E_40, 'o-');
% ylabel('$L_2$ error');
% xlabel('$Refinement$');
% legend('20x40', 'Location', 'best');

% figure;
% loglog(ref, E_40, 'o-');
% ylabel('$L_2$ error');
% xlabel('$Refinement$');
% legend('40x80', 'Location', 'best');

% figure;
% loglog(ref, E_60, 'o-');
% ylabel('$L_2$ error');
% xlabel('$Refinement$');
% legend('60x120', 'Location', 'best');

% pause;
% close all;
% end

% figure;
% semilogy(steps, E_1(:, 1), 'o-'); hold on;
% semilogy(steps, E_1(:, 2), 'o-');
% semilogy(steps, E_old(:, 1), 'o-'); hold on;
% semilogy(steps, E_old(:, 2), 'o-');
% ylabel('$L_2$ error');
% legend('20x40', '40x80', '20x40 old', '40x80 old', 'Location', 'best');
% title('Old BC');

% E_1 = mean(E_1);
% E_2 = mean(E_2);
% E_3 = mean(E_3);
% E_4 = mean(E_4);
% E_old = mean(E_old);

% % do linear fit of E_1 and E_old
% res = [20 40];
% p1 = polyfit(log(res), log(E_1), 1);
% p2 = polyfit(log(res), log(E_old), 1);

% figure;
% loglog(res, E_1, 'o-'); hold on;
% loglog(res, E_old, 'o-');
% ylabel('$L_2$ error');
% xlabel('$N_{particles}$');
% legend('New BC', 'Old BC', 'Location', 'best');

% figure;
% loglog(res, E_1, 'o-'); hold on;
% loglog(res, exp(polyval(p1, log(res))), 'k--');
% loglog(res, E_old, 'o-');
% loglog(res, exp(polyval(p2, log(res))), 'r--');
% text(res(2) * 0.9, E_1(2), ['Slope: ' num2str(p1(1)) '         '], 'HorizontalAlignment', 'right');
% text(res(2) * 1.1, E_old(2), ['         Slope: ' num2str(p2(1))], 'HorizontalAlignment', 'left');
% ylabel('$L_2$ error');
% xlabel('$N_{particles}$');
% legend('New BC', 'New BC linear fit', 'Old BC', 'Old BC linear fit', 'Location', 'best');

function [ux, uy] = Poiseuille_Analytical(y, params)
    g = params(1);
    nu = params(2);
    L = params(3);
    ux = (g / (2 * nu)) .* y .* (L - y);
    uy = zeros(size(y));
end

function u = prof_a_pouiseuille(x, t)

    g = 0.1;
    nu = 0.1;
    L = 1;
    Nterms = 20;
    u = (g / (2 * nu)) .* x .* (L - x);

    for n = 0:Nterms
        u = u - ((4 * g * L ^ 2) / (nu * pi ^ 3 * (2 * n + 1) ^ 3)) .* sin(pi * x * (2 * n + 1) / L) .* exp(- ((2 * n + 1) ^ 2 * pi ^ 2 * nu * t) / (L ^ 2));
    end

end

function [L2] = P_Error(Velocity, ux, uy)
    NormSPH = sqrt(Velocity(:, 1) .^ 2 + Velocity(:, 2) .^ 2);
    NormAnalytical = sqrt(ux .^ 2 + uy .^ 2);

    L2 = sqrt((1 / length(Velocity(:, 1))) * sum ((NormAnalytical - NormSPH) .^ 2));
end

function [L2] = X_Error(Velocity, ux, uy)
    NormSPH = sqrt(Velocity(:, 1) .^ 2);
    NormAnalytical = sqrt(ux .^ 2);

    L2 = sqrt((1 / length(Velocity(:, 1))) * sum ((NormAnalytical - NormSPH) .^ 2));
end

function [L2] = Y_Error(Velocity, ux, uy)
    NormSPH = sqrt(Velocity(:, 2) .^ 2);
    NormAnalytical = sqrt(uy .^ 2);

    L2 = sqrt((1 / length(Velocity(:, 1))) * sum ((NormAnalytical - NormSPH) .^ 2));
end

function [errors] = ComputeErrorAndPlotPoiseuille(datasets, plot_setting)

    nres = length(datasets);
    errors = zeros(nres, 2);

    g = 0.1;
    nu = 0.1;
    L = 1.0;
    params = [g nu L];

    if (plot_setting == 1)
        fig1 = figure;
        % fig2 = figure;
    end

    for k = 1:nres
        data = datasets{k};

        pos = data.Position;
        vel = data.Velocity;
        velT = data.VelocityTransport;

        % Remove boundary particles

        pos = pos(data.Type == 1, :);
        vel = vel(data.Type == 1, :);
        velT = velT(data.Type == 1, :);

        [ux, uy] = Poiseuille_Analytical(pos(:, 2), params);
        errors(k, 1) = P_Error(vel, ux, uy);
        errors(k, 2) = P_Error(velT, ux, uy);

        % fine grid for analytical solution plot
        yfine = linspace(0, 1, 1000);
        [uxfine, uyfine] = Poiseuille_Analytical(yfine, params);

        if (plot_setting == 0)
            continue;
        end

        figure(fig1);
        subplot(2, nres, k);
        hold on;
        plot(yfine, uxfine);
        plot(pos(:, 2), vel(:, 1), '.');
        xlabel('$y$');
        ylabel('$u_x(y)$');
        title([data.PlotName]);
        % axis();
        subplot(2, nres, k + nres);
        hold on;
        plot(yfine, uyfine);
        plot(pos(:, 2), vel(:, 2), '.');
        ylabel('$u_y(y)$');
        xlabel('$y$');
        % axis([1.5 2 -0.1 0.1]);
        title([" $L_2$ error: " + num2str(errors(k, 1))]); %"$L_2$ error (transport): " + num2str(errors(k, 2))])

        % figure(fig2);
        % subplot(2, nres, k);
        % hold on;
        % plot(pos(:, 2), abs(vel(:, 1) - ux), '.');
        % xlabel('$y$');
        % ylabel('$E_x(y)$');
        % % axis([1.5 2 0 0.1]);
        % title([data.PlotName]);
        % subplot(2, nres, k + nres);
        % hold on;
        % plot(pos(:, 2), abs(vel(:, 2) - uy), '.');
        % xlabel('$y$');
        % ylabel('$E_y(y)$');
        % % axis([1.5 2 0 0.1]);
        % title([" $L_2$ error: " + num2str(errors(k, 1)) "$L_2$ error (transport): " + num2str(errors(k, 2))])

    end

    errors = errors(:, 1);
end
