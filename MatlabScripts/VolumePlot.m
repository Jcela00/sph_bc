clear all; close all; clc;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);

% This is for exportgraphics
set(groot, 'defaultFigureUnits', 'centimeters');
% set(groot, 'defaultFigurePosition', [0, 0, 8.5, 6.0]); %single column
set(groot, 'defaultFigurePosition', [0, 0, 16.0, 10.0]); %double column

% MarkerParticleSize = 12;
% dx = 0.1;
% kappa_max = 1 / dx;
% kappa_max_2 = 1 / (2 * dx);
% kappa_max_3 = 1 / (3 * dx);
% kappa_grid = linspace(0, 2 * kappa_max, 1000);
% vol = curvedVol(kappa_grid, dx);
% vref = dx ^ 2;

% figure; hold on;
% plot(kappa_grid, vol(1, :), 'Color', [0.6350 0.0780 0.1840]);
% plot(kappa_grid, vol(2, :), 'Color', [0.4660 0.6740 0.1880]);
% plot(kappa_grid, vol(3, :), 'Color', [0 0.4470 0.7410]);
% xline(kappa_max_3, '--');
% xline(kappa_max_2, '--');
% xline(kappa_max, '--');
% yline(vref, '--');
% xlabel('$\kappa$');
% ylabel('$V$');
% legend('first', 'second', 'third', 'Location', 'best', 'box', 'off');

% xticks([0, kappa_max_3, kappa_max_2, kappa_max]);
% xticklabels({'0', '$\frac{1}{3\Delta x}$', '$\frac{1}{2\Delta x}$', '$\frac{1}{\Delta x}$'});
% yticks([0.01]);
% yticklabels({'$\Delta x^2$'});

% set(gca, 'FontSize', 11); % Adjust axes font size
% set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
% exportgraphics(gcf, 'LatexFigures/VolumeKappa.pdf', 'ContentType', 'vector', 'Resolution', 300);

% r_grid = linspace(dx * 0.1, 10 * dx, 1000);
% inv_r = 1 ./ r_grid;
% vol = curvedVol(inv_r, dx);

% figure; hold on;
% plot(r_grid, vol(1, :), 'Color', [0.6350 0.0780 0.1840]);
% plot(r_grid, vol(2, :), 'Color', [0.4660 0.6740 0.1880]);
% plot(r_grid, vol(3, :), 'Color', [0 0.4470 0.7410]);
% xline(dx, '--');
% xline(2 * dx, '--');
% xline(3 * dx, '--');
% yline(vref, '--');
% xlabel('$R_0$');
% ylabel('$V$');
% legend('first', 'second', 'third', 'Location', 'best', 'box', 'off');

% xticks([dx, 2 * dx, 3 * dx]);
% xticklabels({'$\Delta x$', '$2\Delta x$', '$3\Delta x$'});
% yticks([0.01]);
% yticklabels({'$\Delta x^2$'});

% set(gca, 'FontSize', 11); % Adjust axes font size
% set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
% exportgraphics(gcf, 'LatexFigures/VolumeRadius.pdf', 'ContentType', 'vector', 'Resolution', 300);

MarkerParticleSize = 12;
dx = 0.1;
kappa_max = 1 / dx;
kappa_max_2 = 1 / (2 * dx);
kappa_max_3 = 1 / (3 * dx);
kappa_grid = linspace(-kappa_max, 2 * kappa_max, 1000);
vol = curvedVol(kappa_grid, dx);
vref = dx ^ 2;

figure; hold on;
plot(kappa_grid, vol(1, :), 'Color', [0.6350 0.0780 0.1840]);
plot(kappa_grid, vol(2, :), 'Color', [0.4660 0.6740 0.1880]);
plot(kappa_grid, vol(3, :), 'Color', [0 0.4470 0.7410]);
xline(-kappa_max, '--');
xline(kappa_max_3, '--');
xline(kappa_max_2, '--');
xline(kappa_max, '--');
yline(vref, '--');
xlabel('$\kappa$');
ylabel('$V$');
legend('first', 'second', 'third', 'Location', 'best', 'box', 'off');

xticks([-kappa_max, 0, kappa_max_3, kappa_max_2, kappa_max]);
xticklabels({'$-\frac{1}{\Delta x}$', '0', '$\frac{1}{3\Delta x}$', '$\frac{1}{2\Delta x}$', '$\frac{1}{\Delta x}$'});
yticks([dx ^ 2]);
yticklabels({'$\Delta x^2$'});

set(gca, 'FontSize', 11); % Adjust axes font size
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
exportgraphics(gcf, 'LatexFigures/VolumeKappa.pdf', 'ContentType', 'vector', 'Resolution', 300);

r_grid1 = [linspace(-5 * dx, -0.01 * dx, 1000)];
r_grid2 = [linspace(0.01 * dx, 5 * dx, 1000)];
inv_r1 = 1 ./ r_grid1;
inv_r2 = 1 ./ r_grid2;
vol1 = curvedVol(inv_r1, dx);
vol2 = curvedVol(inv_r2, dx);

figure; hold on;
plot(r_grid1, vol1(1, :), 'Color', [0.6350 0.0780 0.1840]);
plot(r_grid1, vol1(2, :), 'Color', [0.4660 0.6740 0.1880]);
plot(r_grid1, vol1(3, :), 'Color', [0 0.4470 0.7410]);
plot(r_grid2, vol2(1, :), 'Color', [0.6350 0.0780 0.1840]);
plot(r_grid2, vol2(2, :), 'Color', [0.4660 0.6740 0.1880]);
plot(r_grid2, vol2(3, :), 'Color', [0 0.4470 0.7410]);

xline(-dx, '--');
xline(-2 * dx, '--');
xline(-3 * dx, '--');
xline(dx, '--');
xline(2 * dx, '--');
xline(3 * dx, '--');
yline(vref, '--');
xlabel('$R_0$');
ylabel('$V$');
legend('first', 'second', 'third', 'Location', 'best', 'box', 'off');
axis([-inf inf 0 4 * dx ^ 2])
xticks([-3 * dx, -2 * dx, -dx, dx, 2 * dx, 3 * dx]);
xticklabels({'$-3\Delta x$', '$-2\Delta x$', '$-\Delta x$', '$\Delta x$', '$2\Delta x$', '$3\Delta x$'});
yticks([dx ^ 2]);
yticklabels({'$\Delta x^2$'});

set(gca, 'FontSize', 11); % Adjust axes font size
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11); % Apply to all text in the figure
exportgraphics(gcf, 'LatexFigures/VolumeRadius.pdf', 'ContentType', 'vector', 'Resolution', 300);

function [vol] = curvedVol(kappa, dx)
    ncol = size(kappa, 2);
    vol = zeros(3, ncol);

    for n = 1:3
        vol(n, :) = 0.5 * (2 * dx + dx * dx * kappa - 2 * n * dx * dx * kappa) * dx;
    end

    vol = max(0, vol);

end

function R = Volume2Radius(Volume)
    R = sqrt(Volume / pi);
end
