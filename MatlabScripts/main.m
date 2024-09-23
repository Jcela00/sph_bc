clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 1440 1440]);

% sph parameters
rho0 = 1;
dim = 2;
filenum_end = 200;
Nfluid = [60 40 1];
Nboundary = [1, 0, 0];
Hconst = 1.0;

% 0 for triangle
% 1 for cylinder
scenario = 0;

if scenario == 0
    dp = 0.025;
    triangle = ParticleData('../CSV_Data/tri', 0, ['triangle'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    x = linspace(0.3, 1.0, 100);
    y = linspace(0.3, 0.8, 100);

elseif scenario == 1
    dp = 0.002;
    triangle = ParticleData('../CSV_Data/cyl', 0, ['triangle'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
    x = linspace(0.025, 0.075, 100);
    y = linspace(0.025, 0.075, 100);
end

r = triangle.Position(triangle.Type == 0, :);
normal = triangle.Normal(triangle.Type == 0, :);

if scenario == 0
    % remove particles with r(2) < dp or r(2) > 1 - dp

    normal = normal(r(:, 2) > dp & r(:, 2) < 1 - dp, :);
    r = r(r(:, 2) > dp & r(:, 2) < 1 - dp, :);

end

triangle.Position = r;
triangle.Normal = normal;
triangle.Type = zeros(length(r(:, 1)), 1);
r = r(:, 1:2);
normal = normal(:, 1:2);

[r, idx] = sortrows(r);
normal = normal(idx, :);

for n = 1:length(r(:, 1))

    fig1 = figure; hold on;

    for i = 1:length(x)
        xx = x(i);

        for j = 1:length(y)
            yy = y(j);

            if (~isInside([xx yy], scenario))
                color = 'k';

                % [bool, dist2marker, bool1, bool2, bool3] = checkCriterion_m_3(r(n, :), normal(n, :), [xx yy], dp);
                [bool, dist2marker, bool1, bool2, bool3] = checkCriterion_m_3(r(n, :), normal(n, :), [xx yy], dp);

                if bool == true
                    color = 'r';
                end

                if (bool && bool1)
                    color = 'g';
                end

                if (bool && bool2)
                    color = 'b';
                end

                if (bool && bool3)
                    color = 'm';
                end

                plot(xx, yy, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerSize', 2);
            end

        end

    end

    triangle.plotNormals(fig1);
    % plot(r(n, 1), r(n, 2), 'o', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'MarkerSize', 10);
    plotCircle(fig1, r(n, :), 3.0 * dp, 'y');
    plot([r(n, 1) r(n, 1) - 2.5 * dp * normal(n, 1)], [r(n, 2) r(n, 2) - 2.5 * dp * normal(n, 2)], '--k', 'LineWidth', 2);
    plot(r(n, 1) - 0.5 * dp * normal(n, 1), r(n, 2) - 0.5 * dp * normal(n, 2), 'o', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
    plotCircle(fig1, r(n, :) - 0.5 * dp * normal(n, :), 3.0 * dp, 'g');
    plot(r(n, 1) - 1.5 * dp * normal(n, 1), r(n, 2) - 1.5 * dp * normal(n, 2), 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    plotCircle(fig1, r(n, :) - 1.5 * dp * normal(n, :), 3.0 * dp, 'b');
    plot(r(n, 1) - 2.5 * dp * normal(n, 1), r(n, 2) - 2.5 * dp * normal(n, 2), 'o', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'MarkerSize', 6);
    plotCircle(fig1, r(n, :) - 2.5 * dp * normal(n, :), 3.0 * dp, 'm');
    axis equal;
    % save as png
    % saveas(fig1, ['figs/Criterion_full_', num2str(n), '.png']);
    % save as pdf
    exportgraphics(gcf, ['figs/checkCriterion_m3_', num2str(n), '.pdf'], 'ContentType', 'vector');

    close all;
end

function [bool, dist2marker, bool1, bool2, bool3] = checkCriterion_nothing(Q, normal, P, dp)

    % Q position of marker particle with normal normal, P position of the fluid particle
    dist2marker = distance_own(Q, P);
    dist1 = distance_own(Q - 0.5 * normal * dp, P);
    dist2 = distance_own(Q - 1.5 * normal * dp, P);
    dist3 = distance_own(Q - 2.5 * normal * dp, P);

    if dist2marker < 3.0 * dp
        bool = true;
    else
        bool = false;
    end

    if dist1 < 3.0 * dp
        bool1 = true;
    else
        bool1 = false;
    end

    if dist2 < 3.0 * dp
        bool2 = true;
    else
        bool2 = false;
    end

    if dist3 < 3.0 * dp
        bool3 = true;
    else
        bool3 = false;
    end

end

function [bool, dist2marker, bool1, bool2, bool3] = checkCriterion_m_3(Q, normal, P, dp)

    % Q position of marker particle with normal normal, P position of the fluid particle
    dist2marker = distance_own(Q, P);
    dist1 = distance_own(Q - 0.5 * normal * dp, P);
    dist2 = distance_own(Q - 1.5 * normal * dp, P);
    dist3 = distance_own(Q - 2.5 * normal * dp, P);

    bool = dist2marker < dist3;
    bool1 = dist1 < 3.0 * dp;
    bool2 = dist2 < 3.0 * dp;
    bool3 = dist3 < 3.0 * dp;

end

function [bool, dist2marker, bool1, bool2, bool3] = checkCriterion_lwall(Q, normal, P, dp)

    % Q position of marker particle with normal normal, P position of the fluid particle
    dist2marker = distance_own(Q, P);

    lwall = [0.5 * dp 1.5 * dp 2.5 * dp];

    if dist2marker < 3.0 * dp
        bool = true;
    else
        bool = false;
    end

    bool1 = dist2marker + lwall(1) < 3.0 * dp;
    bool2 = dist2marker + lwall(2) < 3.0 * dp;
    bool3 = dist2marker + lwall(3) < 3.0 * dp;

end

function [bool, dist2marker, bool1, bool2, bool3] = checkCriterion_lwallmod(Q, normal, P, dp)

    % Q position of marker particle with normal normal, P position of the fluid particle
    dist2marker = distance_own(Q, P);

    r0 = sqrt(3.0 ^ 2 - (0.5 + 0) ^ 2) * dp;
    r1 = sqrt(3.0 ^ 2 - (0.5 + 1) ^ 2) * dp;
    r2 = sqrt(3.0 ^ 2 - (0.5 + 2) ^ 2) * dp;

    if dist2marker < 3.0 * dp
        bool = true;
    else
        bool = false;
    end

    bool1 = dist2marker < r0;
    bool2 = dist2marker < r1;
    bool3 = dist2marker < r2;

end

function plotCircle(fig, centre, radius, linestyle)

    th = linspace(0, 2 * pi, 20);
    xunit = radius * cos(th) + centre(1);
    yunit = radius * sin(th) + centre(2);
    figure(fig);
    plot(xunit, yunit, linestyle, 'LineWidth', 1);

end

function dist = distance_own(P, Q)
    dist = sqrt((P(1) - Q(1)) ^ 2 + (P(2) - Q(2)) ^ 2);
end

function boolean = isInside(P, scenario)

    if scenario == 0
        A = [0.4366667 0.4];
        B = [0.9066667 0.4];
        C = [0.9066667 0.7];

        %      C
        %    /  |
        %   /   |
        % A-----B

        boolean = false;
        % check if P is inside the triangle ABC

        % first check if P is inside the rectangle formed by AC
        if (P(1) >= A(1) && P(1) <= C(1) && P(2) >= A(2) && P(2) <= C(2))
            % check if P is below the line AC
            if (P(2) <= (C(2) - A(2)) / (C(1) - A(1)) * (P(1) - A(1)) + A(2))
                boolean = true;
            end

        end

    elseif scenario == 1
        R = 0.02;
        C = [0.05 0.05];

        if (sqrt((P(1) - C(1)) ^ 2 + (P(2) - C(2)) ^ 2) <= R)
            boolean = true;
        else
            boolean = false;
        end

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% % Analytical profile
% tend = 80;
% xfine = linspace(0, 1, 1000);
% ufine_p = prof_a_pouiseuille(xfine, tend);
% ufine_c = prof_a_couette(xfine, tend);
% xsamples = 100; % sample points to evaluate SPH sum
% z_coord = 0.1; % channel z point to evaluate profile

% Poiseuille_NEW_BC_Quintic_Differential = ParticleData('../CSV_Data/Poiseuille_NEW_BC_Quintic_Differential_Re1.250000_40_19/file', filenum_end, ['New BC Differential'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
% Poiseuille_NEW_BC_Quintic_Summation = ParticleData('../CSV_Data/Poiseuille_NEW_BC_Quintic_Summation_Re1.250000_40_19/file', filenum_end, ['New BC Summation'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
% Poiseuille_OLD_BC_Quintic_Differential = ParticleData('../CSV_Data/Poiseuille_OLD_BC_Quintic_Differential_Re1.250000_40_19/file', filenum_end, ['Old BC Differential'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
% Poiseuille_OLD_BC_Quintic_Summation = ParticleData('../CSV_Data/Poiseuille_OLD_BC_Quintic_Summation_Re1.250000_40_19/file', filenum_end, ['Old BC Sumation'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
% obj1 = Poiseuille_NEW_BC_Quintic_Differential;
% obj2 = Poiseuille_NEW_BC_Quintic_Summation;
% obj3 = Poiseuille_OLD_BC_Quintic_Differential;
% obj4 = Poiseuille_OLD_BC_Quintic_Summation;

% fig1 = figure; hold on;
% plot(xfine, ufine_p, 'k', 'DisplayName', 'Analytical Poiseuille'); legend;
% obj1.ScatterParticlePlot(fig1, 'rs', 10);
% obj3.ScatterParticlePlot(fig1, 'b+', 6);

% fig2 = figure; hold on;
% plot(xfine, ufine_p, 'k', 'DisplayName', 'Analytical Poiseuille'); legend;
% obj2.ScatterParticlePlot(fig2, 'rs', 10);
% obj4.ScatterParticlePlot(fig2, 'b+', 6);

% error computation
% Err = zeros(obj1.NpartFluid, 4);
% Position_rectified = zeros(obj1.NpartFluid, 4);
% obj_vec = {obj1, obj2, obj3, obj4};
% ErrTotal = zeros(1, 4);

% for j = 1:4
%     obj = obj_vec{j};

%     for k = 1:obj.Npart

%         if (obj.Type(k) == 1)
%             x = obj.Position(k, 1);
%             u = prof_a_pouiseuille(x, tend);
%             Err(k, j) = (obj.Velocity(k, 1) - u) ^ 2;
%             Position_rectified(k, j) = x;
%         end

%         ErrTotal(j) = sqrt(sum(Err(:, j)) / obj.Npart);
%     end

% end

% figure; hold on;
% semilogy(Position_rectified(:, 1), Err(:, 1), 'rs', 'DisplayName', 'New BC Differential');
% semilogy(Position_rectified(:, 2), Err(:, 2), 'bo', 'DisplayName', 'New BC Summation');
% semilogy(Position_rectified(:, 3), Err(:, 3), 'g+', 'DisplayName', 'Old BC Differential');
% semilogy(Position_rectified(:, 4), Err(:, 4), 'kx', 'DisplayName', 'Old BC Summation');

% legend

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytical velocity profiles

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

function u = prof_a_couette(x, t)

    L = 1;
    V0 = 1.25;
    nu = 0.1;
    Nterms = 10;
    u = V0 * x / L;

    for n = 1:Nterms
        u = u + ((2 * V0 / (n * pi)) * (-1) ^ n) .* sin(n * pi * x / L) .* exp(-nu * n ^ 2 * pi ^ 2 * t / (L ^ 2));
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
