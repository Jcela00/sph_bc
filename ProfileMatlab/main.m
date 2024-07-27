clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 1440 1440]);

% sph parameters
dp = 0.025;
rho0 = 1;
dim = 2;
filenum_end = 200;
Nfluid = [60 40 1];
Nboundary = [1, 0, 0];
Hconst = 1.0;

%proc1 = ParticleData('../CSV_Data/1prc/file', 0, ['1proc'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
%proc2 = ParticleData('../CSV_Data/2prc/file', 0, ['2proc'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
% triangle = ParticleData('../CSV_Data/Triangle/file', 0, ['triangle'], dp, Nfluid, Nboundary, rho0, dim, Hconst);
triangle = ParticleData('trianglesym', 0, ['triangle'], dp, Nfluid, Nboundary, rho0, dim, Hconst);

fig1 = figure;

triangle.plotNormals(fig1); axis equal;

% for k = 0:filenum_end

%     proc1 = proc1.Update(k);
%     proc2 = proc2.Update(k);
%     [pos1, I1] = sortrows(proc1.Position);
%     [pos4, I2] = sortrows(proc2.Position);
%     proc1 = proc1.ReorderData(I1);
%     proc2 = proc2.ReorderData(I2);

%     [objDiff, pack] = CompDiff(proc1, proc2);

%     PosDiff = abs(proc1.Position - proc2.Position);

%     error_rho = pack(1);
%     error_p = pack(2);
%     error_Force = pack(3);
%     error_Velocity = pack(4);
%     error_ForceTransport = pack(5);
%     error_VelocityTransport = pack(6);
%     error_Normal = pack(7);
%     error_Curvature = pack(8);
%     error_ArcLength = pack(9);
%     error_InteractCount = pack(10);
%     number_error = pack(11);

%     if (error_rho ~= 0)
%         fig1 = figure; hold on;
%         objDiff.plotMagnitude(objDiff.Density, '$Density$ ', fig1, k)
%     end

%     if (error_p ~= 0)
%         fig2 = figure; hold on;
%         objDiff.plotMagnitude(objDiff.Pressure, '$Pressure$ ', fig2, k)
%     end

%     if (error_Force ~= 0)
%         fig3 = figure; hold on;
%         objDiff.plotMagnitude(objDiff.Force, '$Force$ ', fig3, k)
%     end

%     if (error_Velocity ~= 0)
%         fig4 = figure; hold on;
%         objDiff.plotMagnitude(objDiff.Velocity, '$Velocity$ ', fig4, k)
%     end

%     if (error_ForceTransport ~= 0)
%         fig5 = figure; hold on;
%         objDiff.plotMagnitude(objDiff.ForceTransport, '$Force Transport$ ', fig5, k)
%     end

%     if (error_VelocityTransport ~= 0)
%         fig6 = figure; hold on;
%         objDiff.plotMagnitude(objDiff.VelocityTransport, '$Velocity Transport$ ', fig6, k)
%     end

%     if (error_Normal ~= 0)
%         fig7 = figure; hold on;
%         objDiff.plotMagnitude(objDiff.Normal, '$Normal$ ', fig7, k)
%     end

%     if (error_Curvature ~= 0)
%         fig8 = figure; hold on;
%         objDiff.plotMagnitude(objDiff.Curvature, '$Curvature$ ', fig8, k)
%     end

%     if error_InteractCount ~= 0
%         fig9 = figure; hold on;
%         objDiff.plotMagnitude(objDiff.InteractCount, '$Interact Count$ ', fig9, k)
%     end

%     if (error_ArcLength ~= 0)
%         fig9 = figure; hold on;
%         objDiff.plotMagnitude(objDiff.ArcLength, '$Arc Length$ ', fig9, k)
%     end

%     disp(['Time Step: ', num2str(k)]);
%     % pack = [totalE_rho, totalE_p, totalE_Force, totalE_Velocity, totalE_ForceTransport, totalE_VelocityTransport, totalE_Normal, totalE_Curvature, totalE_ArcLength];
%     disp(['Max Position Difference: ', num2str(max(norm(PosDiff)))]);
%     disp(['Total Error in Density: ', num2str(error_rho)]);
%     disp(['Total Error in Pressure: ', num2str(error_p)]);
%     disp(['Total Error in Force: ', num2str(error_Force)]);
%     disp(['Total Error in Velocity: ', num2str(error_Velocity)]);
%     disp(['Total Error in Force Transport: ', num2str(error_ForceTransport)]); ;
%     disp(['Total Error in Velocity Transport: ', num2str(error_VelocityTransport)]);
%     disp(['Total Error in Normal: ', num2str(error_Normal)]);
%     disp(['Total Error in Curvature: ', num2str(error_Curvature)]);
%     disp(['Total Error in Interact Count: ', num2str(error_InteractCount)]);
%     % disp(['Total Error in Arc Length: ', num2str(error_Arclength)]);
%     disp(['Difference in number of particles: ', num2str(number_error)]);
%     disp(' ');

%     pause()

%     % figure(fig2); clf;

% end
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
