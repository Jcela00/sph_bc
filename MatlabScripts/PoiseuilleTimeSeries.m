clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1);
set(groot, 'defaultLineMarkerSize', 7);
% This is for exportgraphics
set(groot, 'defaultFigureUnits', 'centimeters');
% set(groot, 'defaultFigurePosition', [0, 0, 8.5, 6.0]); %single column
set(groot, 'defaultFigurePosition', [0, 0, 16.0, 10.0]); %double column

mrksz = 4;

% sph parameters
rho0 = 1; dim = 2; Hconst = 1;

% Poiseuille flow parameters
g = 0.1; nu = 0.1; L = 1;
params = [g, nu, L];

dirname = 'PoiseuilleCouetteTime/PoiseuilleTimeData';

yfine = linspace(0, 1, 1000);

filenumber = [20 100 200 1000];
timesteps = [0.2 1 2 10];

fig2 = figure; hold on;

for k = 1:length(filenumber)

    nfile = filenumber(k);
    data30 = ParticleData(['../CSV_Data/' dirname '/Poiseuille_new_differential_12_30_1rf_1prc/file'], nfile, ['30'], rho0, dim, Hconst);
    data90 = ParticleData(['../CSV_Data/' dirname '/Poiseuille_new_differential_36_90_1rf_1prc/file'], nfile, ['90'], rho0, dim, Hconst);

    [time30, ycoord30, vel30] = extractData(data30);
    [time90, ycoord90, vel90] = extractData(data90);

    ycoord90 = ycoord90(2:3:end);
    vel90 = vel90(2:3:end);

    ua = prof_a_pouiseuille(yfine, time30, params);
    txt = ['t = ' num2str(timesteps(k))];
    text(0.5, max(ua) * 1.01, txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

    plot(yfine, ua, 'k', 'HandleVisibility', 'off');

    if k == length(filenumber)
        figure(fig2);
        plot(ycoord90, vel90, 'ob', 'DisplayName', 'SPH $N_y = 90$');
        plot(ycoord30, vel30, 'xk', 'DisplayName', 'SPH $N_y = 30$');
    else
        figure(fig2);
        plot(ycoord90, vel90, 'ob', 'HandleVisibility', 'off');
        plot(ycoord30, vel30, 'xk', 'HandleVisibility', 'off');
    end

    if (k == length(filenumber))
        fig1 = figure; hold on;
        data30.PlotParticles(fig1, mrksz);
        axis equal; axis tight;
        xlabel('$x$'); ylabel('$y$');
        set(gca, 'FontSize', 11); % Adjust axes font size
        set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
        exportgraphics(gcf, 'LatexFigures/PoiseuilleParticles.pdf', 'ContentType', 'vector', 'Resolution', 300);
    end

end

figure(fig2);
xlabel('$y$'); ylabel('$u_x$');
legend('Location', 'northwest', 'box', 'off');
axis([0 1 0 0.14]);
set(gca, 'FontSize', 11); % Adjust axes font size
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
exportgraphics(gcf, 'LatexFigures/Poiseuille.pdf', 'ContentType', 'vector', 'Resolution', 300);

% Couette flow parameters
v0 = 0.125; nu = 0.1; L = 1;
params = [v0, nu, L];

dirname = 'PoiseuilleCouetteTime/CouetteTimeData';

yfine = linspace(0, 1, 1000);

filenumber = [20 100 200 1000];
timesteps = [0.2 1 2 10];

txtx = [0.772 0.633 0.569 0.528];
txty = [0.03178 0.05138 0.06034 0.06608];

fig3 = figure; hold on;

for k = 1:length(filenumber)

    nfile = filenumber(k);
    data30 = ParticleData(['../CSV_Data/' dirname '/Couette_new_differential_12_30_1rf_1prc/file'], nfile, ['30'], rho0, dim, Hconst);
    data90 = ParticleData(['../CSV_Data/' dirname '/Couette_new_differential_36_90_1rf_1prc/file'], nfile, ['90'], rho0, dim, Hconst);

    [time30, ycoord30, vel30] = extractData(data30);
    [time90, ycoord90, vel90] = extractData(data90);

    % take every third point in the 90 particles, starting from the second
    ycoord90 = ycoord90(2:3:end);
    vel90 = vel90(2:3:end);

    ua = prof_a_couette(yfine, time30, params);
    txt = ['t = ' num2str(timesteps(k))];

    plot(yfine, ua, 'k', 'HandleVisibility', 'off');

    if k == length(filenumber)
        figure(fig3);
        plot(ycoord90, vel90, 'ob', 'DisplayName', 'SPH $N_y = 90$');
        plot(ycoord30, vel30, 'xk', 'DisplayName', 'SPH $N_y = 30$');
        text(txtx(k), txty(k), txt, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    else
        figure(fig3);
        plot(ycoord90, vel90, 'ob', 'HandleVisibility', 'off');
        plot(ycoord30, vel30, 'xk', 'HandleVisibility', 'off');
        text(txtx(k), txty(k), txt, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    end

    if (k == length(filenumber))
        fig4 = figure; hold on;
        data30.PlotParticles(fig4, mrksz);
        axis equal; axis tight;
        xlabel('$x$'); ylabel('$y$');
        set(gca, 'FontSize', 11); % Adjust axes font size
        set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
        exportgraphics(gcf, 'LatexFigures/CouetteParticles.pdf', 'ContentType', 'vector', 'Resolution', 300);
    end

end

% slope = -0.14/1;
% yintercept = 0.14;

% f = @(y) slope * y + yintercept;
% yvals = 0:0.001:1;
% tt = f(yvals);

% plot(yvals, tt, 'r');

figure(fig3);
xlabel('$y$'); ylabel('$u_x$');
lgd = legend('Location', 'northwest', 'box', 'off');
axis([0 1 0 0.14]);
set(gca, 'FontSize', 11); % Adjust axes font size
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
exportgraphics(gcf, 'LatexFigures/Couette.pdf', 'ContentType', 'vector', 'Resolution', 300);

% fig1 = figure; hold on; axis equal;
% data30.PlotParticles(fig1);

% FUNCTIONS
function [time, ycoord, vel] = extractData(ParticleDataSet)
    time = ParticleDataSet.Time;

    pos = ParticleDataSet.Position;
    vel = ParticleDataSet.Velocity;
    velT = ParticleDataSet.VelocityTransport;

    % Remove boundary particles
    pos = pos(ParticleDataSet.Type == ParticleDataSet.TypeFluid, :);
    vel = vel(ParticleDataSet.Type == ParticleDataSet.TypeFluid, :);
    velT = velT(ParticleDataSet.Type == ParticleDataSet.TypeFluid, :);

    % Remove z component (all zeros)
    pos = pos(:, 1:2);
    vel = vel(:, 1:2);
    velT = velT(:, 1:2);

    % sort the particles based on y
    [pos, I] = sortrows(pos, 2);
    vel = vel(I, :);
    velT = velT(I, :);

    ycoord = pos(:, 2);
    vel = vel(:, 1);

    tolerance = 0.001; % Define a tolerance for grouping close y-values

    y_rounded = round(ycoord / tolerance) * tolerance;
    [unique_y, ~, idx] = unique(y_rounded);
    filtered_v = accumarray(idx, vel, [], @mean); % Take average v for each group

    ycoord = unique_y;
    vel = filtered_v;
end

% function [ux, uy] = Poiseuille_Analytical(y, params)
%     g = params(1);
%     nu = params(2);
%     L = params(3);
%     ux = (g / (2 * nu)) .* y .* (L - y);
%     uy = zeros(size(y));
% end

function u = prof_a_couette(x, t, params)
    v0 = params(1);
    nu = params(2);
    L = params(3);

    Nterms = 100;
    u = v0 * x / L;

    for n = 1:Nterms
        u = u + ((2 * v0 / (n * pi)) * (-1) ^ n) .* sin(n * pi * x / L) .* exp(-nu * n ^ 2 * pi ^ 2 * t / (L ^ 2));
    end

end

function u = prof_a_pouiseuille(x, t, params)

    g = params(1);
    nu = params(2);
    L = params(3);

    Nterms = 100;
    u = (g / (2 * nu)) .* x .* (L - x);

    for n = 0:Nterms
        u = u - ((4 * g * L ^ 2) / (nu * pi ^ 3 * (2 * n + 1) ^ 3)) .* sin(pi * x * (2 * n + 1) / L) .* exp(- ((2 * n + 1) ^ 2 * pi ^ 2 * nu * t) / (L ^ 2));
    end

end

function [L1, L1spatial] = L1_error(VelocityNumeric, VelocityAnalytic)
    L1spatial = abs(VelocityNumeric(:, 1) - VelocityAnalytic(:, 1)) + abs(VelocityNumeric(:, 2) - VelocityAnalytic(:, 2));
    N = length(VelocityNumeric(:, 1));
    L1 = (1 / N) * sum(L1spatial);
end
