clear all; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1);
set(groot, 'defaultLineMarkerSize', 6);
% This is for exportgraphics
set(groot, 'defaultFigureUnits', 'centimeters');
% set(groot, 'defaultFigurePosition', [0, 0, 8.5, 6.0]); %single column
% set(groot, 'defaultFigurePosition', [100, 100, 16.0, 4.0]); %double column

set(groot, 'defaultFigurePosition', [100, 100, 16.0, 10.0]); %double column

%% GAUSSIAN FIT OF a(t) FOR MOVING SQUARE
data = csvread('../Motion_Body.dat', 1, 0);
t = data(:, 1);
a = data(:, 2);
v = data(:, 3);
x = data(:, 4);

% initial guesses for parameters
A = max(a);
sigma = std(a);
mu = t(a == A);

initialGuess = [A, mu, sigma]; % [Amplitude, Mean, Standard Deviation]
gaussianModel = fittype('A*exp(-(x - mu)^2 / (2*sigma^2))', 'independent', 'x', 'dependent', 'y');

options = fitoptions('Method', 'NonlinearLeastSquares', ...
    'StartPoint', initialGuess, ...
    'TolFun', 1e-10, ... % Function tolerance
    'TolX', 1e-10, ... % Parameter tolerance
    'MaxIter', 1000000); % Maximum number of iterations

% Associate options with the fit type
gaussianModel = setoptions(gaussianModel, options);
% Fit the model to the data
[fitres, gof, output] = fit(t, a, gaussianModel);

ymodel = fitres(t);
A = fitres.A;
sigma = fitres.sigma;
mu = fitres.mu;

Aown = 2 * sqrt(2);
sigmaown = 0.25 / sqrt(pi);
muown = mu;
yown = A * exp(- (t - mu) .^ 2 / (2 * sigma ^ 2));

error_own = norm(a - yown) / norm(a);
error_fit = norm(a - ymodel) / norm(a);

fprintf('Error own: %f\n', error_own);
fprintf('Error fit: %f\n', error_fit);

% Display the fit parameters
format long g;
disp("A = ");
disp(fitres.A);
disp("mu = ");
disp(fitres.mu);
disp("sigma = ");
disp(fitres.sigma);

figure; hold on;
plot(t, a, 'r');
plot(t, v, 'b');
legend('$a(t) D/U^2$', '$v(t)/U$', 'Location', 'best');
xlabel('$t$');
axis([0 3 -inf inf]);
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/SquareTimeLaw.pdf'], 'ContentType', 'vector', 'Resolution', 300);

% plot problem geometry
figure; hold on;
BigRectanglePos = [0 0];
BigRectangleSize = [10 5];

SquarePos = [1 2];
SquareSize = [1 1];
rectangle('Position', [BigRectanglePos BigRectangleSize], 'EdgeColor', 'r', 'LineWidth', 1);
rectangle('Position', [SquarePos SquareSize], 'EdgeColor', 'k', 'LineWidth', 1);
plot(1.5, 2.5, 'ok', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
axis equal; axis([-0.2 10.2 -0.2 5.2]);
xlabel('$x$'); ylabel('$y$');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/SquareGeometry.pdf'], 'ContentType', 'vector', 'Resolution', 300);

%%%%%%%% SIMULATION PLOTS
set(groot, 'defaultFigurePosition', [100, 100, 16.0, 4.0]); %double column

Drag50 = ReadDragLift('../CSV_Data/MovingObstacle/Drag/Re50.csv', 1, 1, 0.02);
Drag100 = ReadDragLift('../CSV_Data/MovingObstacle/Drag/Re100.csv', 1, 1, 0.01);
Drag150 = ReadDragLift('../CSV_Data/MovingObstacle/Drag/Re150.csv', 1, 1, 0.00666667);

[DragRef50] = ReadRefDrag('../CSV_Data/MovingObstacle/Drag/Force_Re050.csv');
[DragRef100] = ReadRefDrag('../CSV_Data/MovingObstacle/Drag/Force_Re100.csv');
[DragRef150] = ReadRefDrag('../CSV_Data/MovingObstacle/Drag/Force_Re150.csv');

figure; hold on;
plot(Drag50{1}, Drag50{3}, 'k-', 'DisplayName', 'SPH');
plot(DragRef50{1}, DragRef50{2}, 'r--', 'DisplayName', 'DNS');
legend('Location', 'best', 'box', 'off');
xlabel('$t$'); ylabel('$C_D$');
axis([0 8 -inf inf]);
text(0.5, 0.6, '$Re = 50$', 'Units', 'normalized', 'FontSize', 11, 'HorizontalAlignment', 'center');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/SquareDrag50.pdf'], 'ContentType', 'vector', 'Resolution', 300);

figure; hold on;
plot(Drag100{1}, Drag100{3}, 'k-', 'DisplayName', 'SPH');
plot(DragRef100{1}, DragRef100{2}, 'r--', 'DisplayName', 'DNS');
% legend('Location', 'best','box', 'off');
xlabel('$t$'); ylabel('$C_D$');
axis([0 8 -inf inf]);
text(0.5, 0.6, '$Re = 100$', 'Units', 'normalized', 'FontSize', 11, 'HorizontalAlignment', 'center');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/SquareDrag100.pdf'], 'ContentType', 'vector', 'Resolution', 300);

figure; hold on;
plot(Drag150{1}, Drag150{3}, 'k-', 'DisplayName', 'SPH');
plot(DragRef150{1}, DragRef150{2}, 'r--', 'DisplayName', 'DNS');
xlabel('$t$'); ylabel('$C_D$');
% legend('Location', 'best');
axis([0 8 -inf inf]);
text(0.5, 0.6, '$Re = 150$', 'Units', 'normalized', 'FontSize', 11, 'HorizontalAlignment', 'center');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/SquareDrag150.pdf'], 'ContentType', 'vector', 'Resolution', 300);

%%%%%%%%%%%%%%%%%%%% THIS IS TO CREATE THE INTERPOLATED DATA WHICH IS VERY TIME CONSUMING %%%%%%%%%%%%%%%%%%%%%
set(groot, 'defaultFigurePosition', [100, 100, 16.0, 10.0]); %double column

Nx = 600;
Ny = 300;
xi = linspace(0, 10, Nx);
yi = linspace(0, 5, Ny);

% [~, zmatrix_cell] = MS50.SPHinterpVelocityWithObstacle(xi, yi, MS50.H);
% uinterp = zmatrix_cell{1}';
% vinterp = zmatrix_cell{2}';
% winterp = zmatrix_cell{3}';
% save('MS50interp80.mat', 'uinterp', 'vinterp', 'winterp');
% clear zmatrix_cell uinterp vinterp winterp;

% [~, zmatrix_cell] = MS100.SPHinterpVelocityWithObstacle(xi, yi, MS100.H);
% uinterp = zmatrix_cell{1}';
% vinterp = zmatrix_cell{2}';
% winterp = zmatrix_cell{3}';
% save('MS100interp80.mat', 'uinterp', 'vinterp', 'winterp');
% % clear zmatrix_cell uinterp vinterp winterp;

% [~, zmatrix_cell] = MS150.SPHinterpVelocityWithObstacle(xi, yi, MS150.H);
% uinterp = zmatrix_cell{1}';
% vinterp = zmatrix_cell{2}';
% winterp = zmatrix_cell{3}';
% save('MS150interp80.mat', 'uinterp', 'vinterp', 'winterp');
% clear zmatrix_cell uinterp vinterp winterp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho0 = 1; dim = 2; Hconst = 1;
MS50_50 = ParticleData('../CSV_Data/MovingObstacle/Re50/file', 50, ['Re50'], rho0, dim, Hconst);
MS100_50 = ParticleData('../CSV_Data/MovingObstacle/Re100/file', 50, ['Re100'], rho0, dim, Hconst);
MS150_50 = ParticleData('../CSV_Data/MovingObstacle/Re150/file', 50, ['Re150'], rho0, dim, Hconst);

MS50_80 = ParticleData('../CSV_Data/MovingObstacle/Re50/file', 80, ['Re50'], rho0, dim, Hconst);
MS100_80 = ParticleData('../CSV_Data/MovingObstacle/Re100/file', 80, ['Re100'], rho0, dim, Hconst);
MS150_80 = ParticleData('../CSV_Data/MovingObstacle/Re150/file', 80, ['Re150'], rho0, dim, Hconst);

Datasets = {MS50_50, MS100_50, MS150_50, MS50_80, MS100_80, MS150_80};
filenames = {'MS50interp50.mat', 'MS100interp50.mat', 'MS150interp50.mat', 'MS50interp80.mat', 'MS100interp80.mat', 'MS150interp80.mat'};
names = {'Re50_50_grid', 'Re100_50_grid', 'Re150_50_grid',
         'Re50_80_grid', 'Re100_80_grid', 'Re150_80_grid'};

names = {'Re50_50', 'Re100_50', 'Re150_50',
         'Re50_80', 'Re100_80', 'Re150_80'};

% Datasets = {MS50_50};
% filenames = {'MS50interp50.mat'};
% names = {'Re50_50'};

for k = 1:length(Datasets)
    PlotMovingSquareContours(Datasets{k}, filenames{k}, names{k});
end

function [] = PlotMovingSquareContours(Dataset, interp_filename, name)

    Nx = 600;
    Ny = 300;
    xi = linspace(0, 10, Nx);
    yi = linspace(0, 5, Ny);

    % dx = xi(2) - xi(1);
    % dy = yi(2) - yi(1);

    % extract reynolds number and time from name
    % name is formated as ReXXX_YY where XXX is the Reynolds number and YY is the number of particles
    splitName = split(name, '_');
    Re = str2double(splitName{1}(3:end));
    N = str2double(splitName{2});

    if (N == 50)
        time = 5;
    elseif (N == 80)
        time = 8;
    end

    txt = {['$Re = ' num2str(Re) '$'], ['$t = ' num2str(time) '$']};

    load(['dotMatData/' interp_filename]);

    % [du_dx, du_dy] = gradient(uinterp, dy, dx);
    % [dv_dx, dv_dy] = gradient(vinterp, dy, dx);

    % winterp = dv_dx - du_dy;

    squareCentre = Dataset.ForceTransport(Dataset.ObstacleIndexes, 1:2);
    squareCentre = squareCentre(1, :);
    squareSide = 1;

    fig3 = figure; hold on;
    levels = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5];
    ncolors = length(levels) + 1;
    cmap = BlueGreenYellowRed(ncolors);
    contourf(xi, yi, sqrt(uinterp .^ 2 + vinterp .^ 2), levels, 'EdgeColor', 'none');
    plotSolidSquare(fig3, squareCentre, squareSide, [0.7 0.7 0.7]);
    colormap(cmap);
    cbar = colorbar;
    set(cbar, 'TickLabelInterpreter', 'latex');
    caxis([0 1.5]);
    axis equal;
    axis([0 10 1 4]);
    xlabel('$x$'); ylabel('$y$');
    text(0.05, 0.5, txt, 'Units', 'normalized', 'FontSize', 11);
    set(gca, 'FontSize', 11);
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
    exportgraphics(gcf, ['LatexFigures/SquareVelMag' name '.pdf'], 'ContentType', 'vector', 'Resolution', 300);

    fig2 = figure; hold on;
    levels = [-inf -6 -5 -4 -3 -2 -1 -0.5 0.5 1 2 3 4 5 6 inf];
    ncolors = length(levels) + 1;
    cmap = BlueYellowRedCmap(ncolors);
    contourf(xi, yi, winterp, levels, 'EdgeColor', 'none');

    % tick_positions = [-4 -3.2 -2.4 -1.6 -0.8 0.8 1.6 2.4 3.2 4];
    % tick_labels = {'-4', '-3.2', '-2.4', '-1.6', '-0.8', '0.8', '1.6', '2.4', '3.2', '4'};

    colormap(cmap);
    cbar = colorbar;
    caxis([-4.9 4.9]);
    % cbar.Ticks = tick_positions;
    % cbar.TickLabels = tick_labels;

    plotSolidSquare(fig2, squareCentre, squareSide, [0.7 0.7 0.7]);
    set(cbar, 'TickLabelInterpreter', 'latex');
    axis equal;
    axis([0 10 1 4]);
    xlabel('$x$'); ylabel('$y$');
    text(0.05, 0.5, txt, 'Units', 'normalized', 'FontSize', 11);
    set(gca, 'FontSize', 11);
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
    exportgraphics(gcf, ['LatexFigures/SquareVorticity' name '.pdf'], 'ContentType', 'vector', 'Resolution', 300);

end

function cmap = BlueYellowRedCmap(nColors)
    % Define adjusted colors
    blue = [0, 0, 1]; % Darker blue
    light_blue = [0, 0.8, 1]; % Shiny cyan blue
    white = [1, 1, 1]; % White
    yellow = [1, 0.8, 0]; % Shiny orange
    red = [1, 0, 0];

    % Split evenly between negative and positive sides
    nPart = ceil(nColors / 4);

    % Interpolate colors
    % nPart = floor(nColors / 4);
    blue_part = [linspace(blue(1), light_blue(1), nPart)', linspace(blue(2), light_blue(2), nPart)', linspace(blue(3), light_blue(3), nPart)';
                 linspace(light_blue(1), white(1), nPart)', linspace(light_blue(2), white(2), nPart)', linspace(light_blue(3), white(3), nPart)'];
    red_part = [linspace(white(1), yellow(1), nPart)', linspace(white(2), yellow(2), nPart)', linspace(white(3), yellow(3), nPart)';
                linspace(yellow(1), red(1), nPart)', linspace(yellow(2), red(2), nPart)', linspace(yellow(3), red(3), nPart)'];

    cmap = [blue_part; red_part];

end

function cmap = BlueGreenYellowRed(nColors)
    white = [1, 1, 1]; % White
    light_blue = [0.6, 0.8, 1]; % Light Blue
    green = [0, 0.8, 0]; % Green
    yellow = [1, 1, 0]; % Yellow
    red = [1, 0, 0]; % Red

    cmap = [linspace(white(1), light_blue(1), ceil(nColors / 4))', linspace(white(2), light_blue(2), ceil(nColors / 4))', linspace(white(3), light_blue(3), ceil(nColors / 4))';
            linspace(light_blue(1), green(1), ceil(nColors / 4))', linspace(light_blue(2), green(2), ceil(nColors / 4))', linspace(light_blue(3), green(3), ceil(nColors / 4))';
            linspace(green(1), yellow(1), ceil(nColors / 4))', linspace(green(2), yellow(2), ceil(nColors / 4))', linspace(green(3), yellow(3), ceil(nColors / 4))';
            linspace(yellow(1), red(1), ceil(nColors / 4))', linspace(yellow(2), red(2), ceil(nColors / 4))', linspace(yellow(3), red(3), ceil(nColors / 4))'];
end

function plotSolidSquare(figHandle, center, D, color)
    % Extract center coordinates
    cx = center(1);
    cy = center(2);

    % Calculate the vertices of the square
    halfD = D / 2;
    x = [cx - halfD, cx + halfD, cx + halfD, cx - halfD];
    y = [cy - halfD, cy - halfD, cy + halfD, cy + halfD];

    figure(figHandle);
    % Plot the square
    fill(x, y, color, 'EdgeColor', 'k'); % Solid fill with black border
end

function [DragDataset] = ReadDragLift(filename, t0, renormalize, nu)
    data = csvread(filename, 1, 0);
    t = data(:, 1) / t0;
    u = data(:, 2);
    drag = data(:, 3);
    lift = data(:, 4);

    if (renormalize == 1)
        drag = -2 * nu * drag .* u;
        lift = -2 * nu * lift .* u;
    end

    DragDataset = {t, u, drag, lift};
end

function [DragDataset] = ReadRefDrag(filename)
    data = csvread(filename, 1, 0);
    t = data(:, 1);
    dragp = data(:, 2);
    dragv = data(:, 3);
    drag = dragp + dragv;

    DragDataset = {t, drag};
end
