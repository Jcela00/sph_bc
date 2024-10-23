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

steps = [1:1:1000];
E_1 = zeros(length(steps), 1);
E_2 = zeros(length(steps), 1);

for k = 1:length(steps)
    nfile = steps(k);
    data1 = ParticleData('../CSV_Data/Poiseuille_new_summation_24_60_1rf_1prc/file', nfile, ['cfl1'], rho0, dim, Hconst);
    data2 = ParticleData('../CSV_Data/Poiseuille_new_summation_unordered_24_60_1rf_1prc/file', nfile, ['unordered'], rho0, dim, Hconst);
    set1 = {data1};
    set2 = {data2};

    plot_setting = 0;

    if (k == length(steps))
        plot_setting = 1;
    end

    E_1(k, 1) = ComputeErrorAndPlotPoiseuille(set1, plot_setting);
    E_2(k, 1) = ComputeErrorAndPlotPoiseuille(set2, plot_setting);
    k
end

figure; hold on;
plot(steps, E_1, 'o-');
plot(steps, E_2, 'o-'); legend('Ordered', 'Unordered');

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
    nu = 0.01;
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

        pos = pos(data.Type == data.TypeFluid, :);
        vel = vel(data.Type == data.TypeFluid, :);
        velT = velT(data.Type == data.TypeFluid, :);

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
