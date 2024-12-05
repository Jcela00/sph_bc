clear all; clc; close all;
mrksz = 8;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1.0);
set(groot, 'defaultLineMarkerSize', mrksz);
set(groot, 'defaultFigureUnits', 'centimeters');
set(groot, 'defaultFigurePosition', [100, 100, 16.0, 16.0]); %double column

rho0 = 1; dim = 2; Hconst = 1;
ConcaveCylinder = ParticleData(['../CSV_Data/CurvaturePlot/ConcaveCylinder/file'], 0, ["ConcaveCylinder"], rho0, dim, Hconst);
ConvexCylinder = ParticleData(['../CSV_Data/CurvaturePlot/ConvexCylinder/file'], 0, ["ConvexCylinder"], rho0, dim, Hconst);
ConcaveSquare = ParticleData(['../CSV_Data/CurvaturePlot/ConcaveSquare/file'], 0, ["ConcaveSquare"], rho0, dim, Hconst);
ConvexSquare = ParticleData(['../CSV_Data/CurvaturePlot/ConvexSquare/file'], 0, ["ConvexSquare"], rho0, dim, Hconst);
Ellipse = ParticleData(['../CSV_Data/CurvaturePlot/FlatEllipse/file'], 0, ["Ellipse"], rho0, dim, Hconst);
Triangle = ParticleData(['../CSV_Data/CurvaturePlot/Triangle/file'], 0, ["Triangle"], rho0, dim, Hconst);
Cavity = ParticleData(['../CSV_Data/CurvaturePlot/Cavity/file'], 0, ["Cavity"], rho0, dim, Hconst);
Cavity.ObstacleIndexes = Cavity.BoundaryIndexes;

datasets = {ConcaveCylinder, ConvexCylinder, ConcaveSquare, ConvexSquare, Ellipse, Triangle, Cavity};
names = ["ConcaveCylinder", "ConvexCylinder", "ConcaveSquare", "ConvexSquare", "Ellipse", "Triangle", "Cavity"];
indexsettings = {0, 0, 1, 0, 0, 0, 0};
vectorsizes = {0.2 0.2 0.2 0.2 0.2 0.2, 0.05};

curvLimitCircle = [-1 1];
curvLimitSquare = [-9 9];
curvLimitEllipse = [0 6];
curvLimitTriangle = [0 10];

curvatureLimits = {curvLimitCircle, curvLimitCircle, curvLimitSquare, curvLimitSquare, curvLimitEllipse, curvLimitTriangle, curvLimitSquare};

datasets = {Cavity};
names = ["Cavity"];
indexsettings = {0};
vectorsizes = {0.04};

curvLimitCircle = [-1 1];
curvLimitSquare = [-9 9];
curvLimitEllipse = [0 6];
curvLimitTriangle = [0 10];

curvatureLimits = {curvLimitSquare};

for k = 1:length(datasets)
    fig = figure; hold on;
    datasets{k}.PlotCurvatureNormals(fig, mrksz, indexsettings{k}, curvatureLimits{k}, vectorsizes{k});
    axis equal;
    exportgraphics(fig, ['LatexFigures/Curvature' names{k} '.pdf'], 'ContentType', 'vector', 'Resolution', 300);

end

% Cylinder45 = ParticleData(['../CSV_Data/CurvaturePlot/Cylinder45/file'], 0, ["Cylinder 45"], rho0, dim, Hconst);
% Cylinder90 = ParticleData(['../CSV_Data/CurvaturePlot/Cylinder90/file'], 0, ["Cylinder 90"], rho0, dim, Hconst);
% Cylinder180 = ParticleData(['../CSV_Data/CurvaturePlot/Cylinder180/file'], 0, ["Cylinder 180"], rho0, dim, Hconst);
% datasets = {Cylinder45, Cylinder90, Cylinder180};

% for k = 1:length(datasets)

%     data = datasets{k};
%     fig1 = figure; hold on;

%     N = length(data.ObstacleIndexes);
%     fprintf('Number of points: %d\n', N);

%     pos = data.Position(data.ObstacleIndexes, 1:2);
%     centre = data.ForceTransport(data.ObstacleIndexes, 1:2);
%     normals = data.Normal(data.ObstacleIndexes, 1:2);
%     % curv = data.Volumes(data.ObstacleIndexes, 2);
%     centre = centre(1, :);
%     error = [];

%     for j = 1:length(pos(:, 1))
%         exact_normal = pos(j, :) - centre;
%         exact_normal = exact_normal / norm(exact_normal);
%         error_x = exact_normal(1) - normals(j, 1);
%         error_y = exact_normal(2) - normals(j, 2);
%         error = [error sqrt(error_x .^ 2 + error_y .^ 2)];
%         % quiver(pos(k, 1), pos(k, 2), exact_normal(1), exact_normal(2), 0.1, 'r', 'LineWidth', 2, 'AutoScale', 'off');
%     end

%     fprintf('Mean error: %f\nMax error: %f\n', mean(error), max(error));
%     figure;
%     plot(error * 100, 'LineWidth', 1.5); title(['Error in normal direction for ' data.BaseName]);
% end

% fig1 = figure; hold on;
% ConvexCylinder.PlotCurvatureNormals(fig1, mrksz, 0, curvLimitCircle);

% pos = ConvexCylinder.Position(ConvexCylinder.ObstacleIndexes, 1:2);
% centre = ConvexCylinder.ForceTransport(ConvexCylinder.ObstacleIndexes, 1:2);
% normals = ConvexCylinder.Normal(ConvexCylinder.ObstacleIndexes, 1:2);
% curv = ConvexCylinder.Volumes(ConvexCylinder.ObstacleIndexes, 2);
% centre = centre(1, :);
% error = [];

% for k = 1:length(pos(:, 1))
%     exact_normal = pos(k, :) - centre; exact_normal = exact_normal / norm(exact_normal);
%     error_x = exact_normal(1) - normals(k, 1);
%     error_y = exact_normal(2) - normals(k, 2);
%     error = [error sqrt(error_x .^ 2 + error_y .^ 2)];
%     quiver(pos(k, 1), pos(k, 2), exact_normal(1), exact_normal(2), 0.1, 'r', 'LineWidth', 2, 'AutoScale', 'off');
% end

% figure;
% plot(error * 100, 'LineWidth', 1.5);
% error = error * 100;
% fprintf('Mean error: %f\nMax error: %f\n', mean(error), max(error));

% fig2 = figure; hold on;
% % Ellipse = ParticleData(['../CSV_Data/CurvaturePlot/Ellipse/file'], 0, ["Ellipse"], rho0, dim, Hconst);

% Ellipse.PlotCurvatureNormals(fig2, mrksz, 0, curvLimitEllipse);
% axis equal;

% pos = Ellipse.Position(Ellipse.ObstacleIndexes, 1:2);
% centre = Ellipse.ForceTransport(Ellipse.ObstacleIndexes, 1:2);
% normals = Ellipse.Normal(Ellipse.ObstacleIndexes, 1:2);
% curv = Ellipse.Volumes(Ellipse.ObstacleIndexes, 2);
% centre = centre(1, :);
% error = [];
% Npoints = length(pos(:, 1))

% kappa = [];

% for k = 1:Npoints
%     a = 1;
%     b = 0.4;
%     exact_normal = [2 * (pos(k, 1) - centre(1)) / a ^ 2, 2 * (pos(k, 2) - centre(2)) / b ^ 2];
%     exact_normal = exact_normal / norm(exact_normal);

%     t = acos((pos(k, 1) - centre(1)) / a);
%     % tt = asin((pos(k, 2) - centre(2)) / b)
%     kappa_exact = a * b / (a ^ 2 * sin(t) ^ 2 + b ^ 2 * cos(t) ^ 2) ^ (3/2);
%     kappa = [kappa kappa_exact];
%     error_x = exact_normal(1) - normals(k, 1);
%     error_y = exact_normal(2) - normals(k, 2);
%     error = [error sqrt(error_x .^ 2 + error_y .^ 2)];
%     quiver(pos(k, 1), pos(k, 2), exact_normal(1), exact_normal(2), 0.1, 'r', 'LineWidth', 2, 'AutoScale', 'off');
% end

% figure;
% plot(kappa); hold on;
% plot(curv); legend('Exact', 'Computed');

% ekappa = abs(kappa' - curv) ./ abs(kappa');

% figure;
% plot(ekappa * 100); title('Error in curvature');

% figure;
% plot(error * 100, 'LineWidth', 1.5);
% error = error * 100;
% fprintf('Mean error: %f\nMax error: %f\n', mean(error), max(error));

% fig2 = figure; hold on;
% Ellipse = ParticleData(['../CSV_Data/CurvaturePlot/EllipseRefined/file'], 0, ["Ellipse"], rho0, dim, Hconst);

% Ellipse.PlotCurvatureNormals(fig2, mrksz, 0, curvLimitEllipse);
% axis equal;

% pos = Ellipse.Position(Ellipse.ObstacleIndexes, 1:2);
% centre = Ellipse.ForceTransport(Ellipse.ObstacleIndexes, 1:2);
% normals = Ellipse.Normal(Ellipse.ObstacleIndexes, 1:2);
% centre = centre(1, :);
% error = [];
% Npoints = length(pos(:, 1))

% for k = 1:Npoints
%     a = 1;
%     b = 0.4;
%     exact_normal = [2 * (pos(k, 1) - centre(1)) / a ^ 2, 2 * (pos(k, 2) - centre(2)) / b ^ 2];
%     exact_normal = exact_normal / norm(exact_normal);

%     error_x = exact_normal(1) - normals(k, 1);
%     error_y = exact_normal(2) - normals(k, 2);
%     error = [error sqrt(error_x .^ 2 + error_y .^ 2)];
%     quiver(pos(k, 1), pos(k, 2), exact_normal(1), exact_normal(2), 0.1, 'r', 'LineWidth', 2, 'AutoScale', 'off');
% end

% figure;
% plot(error * 100, 'LineWidth', 1.5);
% error = error * 100;
% fprintf('Mean error: %f\nMax error: %f\n', mean(error), max(error));
