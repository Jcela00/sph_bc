clear all; close all; clc;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigurePosition', [1440 0 600 600]);

dx = 0.1;
kappa_max = 1 / dx;
kappa_max_2 = 1 / (2 * dx);
kappa_max_3 = 1 / (3 * dx);
kappa_grid = linspace(0, 2 * kappa_max, 1000);
vol = curvedVol(kappa_grid, dx);
vref = dx ^ 2;

% figure; hold on;
% plot(kappa_grid, vol(1, :), 'r-');
% plot(kappa_grid, vol(2, :), 'g-');
% plot(kappa_grid, vol(3, :), 'b-');
% xline(kappa_max_3, '--', 'k=1/(3\Deltax)');
% xline(kappa_max_2, '--', 'k=1/(2\Deltax)');
% xline(kappa_max, '--', 'k=1/\Deltax');
% yline(vref, '--', 'V_{flat}');
% xlabel('$\kappa$');
% ylabel('Volume');
% legend('first', 'second', 'third');

% r_grid = linspace(dx * 0.1, 30 * dx, 1000);
% inv_r = 1 ./ r_grid;
% vol = curvedVol(inv_r, dx);

% figure; hold on;
% plot(r_grid, vol(1, :), 'r-');
% plot(r_grid, vol(2, :), 'g-');
% plot(r_grid, vol(3, :), 'b-');
% xline(dx, '--', '\Deltax');
% xline(2 * dx, '--', '2\Deltax');
% xline(3 * dx, '--', '3\Deltax');
% yline(vref, '--', 'V_{flat}');
% xlabel('$R_0$');
% ylabel('Volume');

R0 = 0.05; dx = 0.01;
kappa = 1 / R0;
vols = curvedVol(kappa, dx);
Rs = Volume2Radius(vols); % radius of the imaginary particles
th = pi / 2;

Npoints = ceil (2 * pi * R0 / dx);
dtheta = 2 * pi / Npoints;

pcentral = R0 * [cos(th), sin(th)];
normal_c = dx * [cos(th), sin(th)];
pcentral_back = pcentral - 3 * normal_c;
virtual_central = [pcentral - 0.5 * normal_c; pcentral - 1.5 * normal_c; pcentral - 2.5 * normal_c];

pright = R0 * [cos(th + dtheta), sin(th + dtheta)];
normal_r = dx * [cos(th + dtheta), sin(th + dtheta)];
pright_back = pright - 3 * normal_r;
virtual_right = [pright - 0.5 * normal_r; pright - 1.5 * normal_r; pright - 2.5 * normal_r];

pleft = R0 * [cos(th - dtheta), sin(th - dtheta)];
normal_l = dx * [cos(th - dtheta), sin(th - dtheta)];
pleft_back = pleft - 3 * normal_l;
virtual_left = [pleft - 0.5 * normal_l; pleft - 1.5 * normal_l; pleft - 2.5 * normal_l];

p_hr = R0 * [cos(th + dtheta / 2), sin(th + dtheta / 2)];
p_hr_back = p_hr - 3 * dx * [cos(th + dtheta / 2), sin(th + dtheta / 2)];

p_hl = R0 * [cos(th - dtheta / 2), sin(th - dtheta / 2)];
p_hl_back = p_hl - 3 * dx * [cos(th - dtheta / 2), sin(th - dtheta / 2)];

line_right = [p_hr(1), p_hr_back(1); p_hr(2), p_hr_back(2)];
line_left = [p_hl(1), p_hl_back(1); p_hl(2), p_hl_back(2)];

p_1hr = R0 * [cos(th + 3 * dtheta / 2), sin(th + 3 * dtheta / 2)];
p_1hr_back = p_1hr - 3 * dx * [cos(th + 3 * dtheta / 2), sin(th + 3 * dtheta / 2)];

p_1hl = R0 * [cos(th - 3 * dtheta / 2), sin(th - 3 * dtheta / 2)];
p_1hl_back = p_1hl - 3 * dx * [cos(th - 3 * dtheta / 2), sin(th - 3 * dtheta / 2)];

line_1right = [p_1hr(1), p_1hr_back(1); p_1hr(2), p_1hr_back(2)];
line_1left = [p_1hl(1), p_1hl_back(1); p_1hl(2), p_1hl_back(2)];

fig1 = figure; hold on;
PlotCircle(fig1, [0, 0], R0, 3, '-');
PlotCircle(fig1, [0, 0], R0 -1 * dx, 1, '--');
PlotCircle(fig1, [0, 0], R0 -2 * dx, 1, '--');
PlotCircle(fig1, [0, 0], R0 -3 * dx, 1, '--');
quiver(pcentral(1), pcentral(2), normal_c(1), normal_c(2), 'k');
quiver(pright(1), pright(2), normal_r(1), normal_r(2), 'k');
quiver(pleft(1), pleft(2), normal_l(1), normal_l(2), 'k');
plot(line_right(1, :), line_right(2, :), ':k');
plot(line_left(1, :), line_left(2, :), ':k');
plot(line_1right(1, :), line_1right(2, :), ':k');
plot(line_1left(1, :), line_1left(2, :), ':k');

scale_factor = 0.5;
PlotCircle(fig1, [virtual_central(1, :)], Rs(1) * scale_factor, 0.5, '-');
PlotCircle(fig1, [virtual_central(2, :)], Rs(2) * scale_factor, 0.5, '-');
PlotCircle(fig1, [virtual_central(3, :)], Rs(3) * scale_factor, 0.5, '-');

PlotCircle(fig1, [virtual_left(1, :)], Rs(1) * scale_factor, 0.5, '-');
PlotCircle(fig1, [virtual_left(2, :)], Rs(2) * scale_factor, 0.5, '-');
PlotCircle(fig1, [virtual_left(3, :)], Rs(3) * scale_factor, 0.5, '-');

PlotCircle(fig1, [virtual_right(1, :)], Rs(1) * scale_factor, 0.5, '-');
PlotCircle(fig1, [virtual_right(2, :)], Rs(2) * scale_factor, 0.5, '-');
PlotCircle(fig1, [virtual_right(3, :)], Rs(3) * scale_factor, 0.5, '-');

vec = [cos(th - 2.0 * dtheta), sin(th - 2.0 * dtheta)];
facc = 1.01;
plabels = [(facc * R0) * vec; facc * (R0 - 1 * dx) * vec; facc * (R0 - 2 * dx) * vec; facc * (R0 - 3 * dx) * vec];
% Add lebels to the volumes
text(virtual_central(1, 1), virtual_central(1, 2), '$V_0$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);
text(virtual_central(2, 1), virtual_central(2, 2), '$V_1$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);
text(virtual_central(3, 1), virtual_central(3, 2), '$V_2$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);

% Add the labels along the curves
text(plabels(1, 1), plabels(1, 2), '$R_0$', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline', 'FontSize', 14);
text(plabels(2, 1), plabels(2, 2), '$R_1$', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline', 'FontSize', 14);
text(plabels(3, 1), plabels(3, 2), '$R_2$', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline', 'FontSize', 14);
text(plabels(4, 1), plabels(4, 2), '$R_3$', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline', 'FontSize', 14);

% Add circular sector for dxwall
PlotCircularSector(fig1, [0, 0], R0, th - dtheta / 2, th + dtheta / 2, 4, '-', 'r');
th_offset = -0.25 * dtheta;
dxlabelpos = 1.02 * R0 * [cos(th + th_offset), sin(th + th_offset)];
text(dxlabelpos(1), dxlabelpos(2), '$\Delta x$', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 14);
% plot marker particles
plot(pcentral(1), pcentral(2), 'ob', 'MarkerFaceColor', 'b');
plot(pright(1), pright(2), 'ob', 'MarkerFaceColor', 'b');
plot(pleft(1), pleft(2), 'ob', 'MarkerFaceColor', 'b');
axis([-2 * dx, 2 * dx, R0 - 4 * dx, R0 + 1 * dx]);
axis equal;
% set(gca, 'Visible', 'off');
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
% print -dpdf -painters epsFig
print(gcf, 'convex.pdf', '-dpdf', '-bestfit');

%%% CONVEX CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
vols = curvedVol(-kappa, dx);
Rs = Volume2Radius(vols); % radius of the imaginary particles
th = 3 * pi / 2;
pcentral = R0 * [cos(th), sin(th)];
normal_c = -dx * [cos(th), sin(th)];
virtual_central = [pcentral - 0.5 * normal_c; pcentral - 1.5 * normal_c; pcentral - 2.5 * normal_c];

pright = R0 * [cos(th + dtheta), sin(th + dtheta)];
normal_r = -dx * [cos(th + dtheta), sin(th + dtheta)];
virtual_right = [pright - 0.5 * normal_r; pright - 1.5 * normal_r; pright - 2.5 * normal_r];

pleft = R0 * [cos(th - dtheta), sin(th - dtheta)];
normal_l = -dx * [cos(th - dtheta), sin(th - dtheta)];
virtual_left = [pleft - 0.5 * normal_l; pleft - 1.5 * normal_l; pleft - 2.5 * normal_l];

p_hr = R0 * [cos(th + dtheta / 2), sin(th + dtheta / 2)];
p_hr_back = p_hr + 3 * dx * [cos(th + dtheta / 2), sin(th + dtheta / 2)];

p_hl = R0 * [cos(th - dtheta / 2), sin(th - dtheta / 2)];
p_hl_back = p_hl + 3 * dx * [cos(th - dtheta / 2), sin(th - dtheta / 2)];

line_right = [p_hr(1), p_hr_back(1); p_hr(2), p_hr_back(2)];
line_left = [p_hl(1), p_hl_back(1); p_hl(2), p_hl_back(2)];

p_1hr = R0 * [cos(th + 3 * dtheta / 2), sin(th + 3 * dtheta / 2)];
p_1hr_back = p_1hr + 3 * dx * [cos(th + 3 * dtheta / 2), sin(th + 3 * dtheta / 2)];

p_1hl = R0 * [cos(th - 3 * dtheta / 2), sin(th - 3 * dtheta / 2)];
p_1hl_back = p_1hl + 3 * dx * [cos(th - 3 * dtheta / 2), sin(th - 3 * dtheta / 2)];

line_1right = [p_1hr(1), p_1hr_back(1); p_1hr(2), p_1hr_back(2)];
line_1left = [p_1hl(1), p_1hl_back(1); p_1hl(2), p_1hl_back(2)];

fig2 = figure; hold on;
PlotCircle(fig2, [0, 0], R0, 3, '-');
PlotCircle(fig2, [0, 0], R0 +1 * dx, 1, '--');
PlotCircle(fig2, [0, 0], R0 +2 * dx, 1, '--');
PlotCircle(fig2, [0, 0], R0 +3 * dx, 1, '--');
quiver(pcentral(1), pcentral(2), normal_c(1), normal_c(2), 'k');
quiver(pright(1), pright(2), normal_r(1), normal_r(2), 'k');
quiver(pleft(1), pleft(2), normal_l(1), normal_l(2), 'k');
plot(line_right(1, :), line_right(2, :), ':k');
plot(line_left(1, :), line_left(2, :), ':k');
plot(line_1right(1, :), line_1right(2, :), ':k');
plot(line_1left(1, :), line_1left(2, :), ':k');

PlotCircle(fig2, [virtual_central(1, :)], Rs(1) * scale_factor, 0.5, '-');
PlotCircle(fig2, [virtual_central(2, :)], Rs(2) * scale_factor, 0.5, '-');
PlotCircle(fig2, [virtual_central(3, :)], Rs(3) * scale_factor, 0.5, '-');

PlotCircle(fig2, [virtual_left(1, :)], Rs(1) * scale_factor, 0.5, '-');
PlotCircle(fig2, [virtual_left(2, :)], Rs(2) * scale_factor, 0.5, '-');
PlotCircle(fig2, [virtual_left(3, :)], Rs(3) * scale_factor, 0.5, '-');

PlotCircle(fig2, [virtual_right(1, :)], Rs(1) * scale_factor, 0.5, '-');
PlotCircle(fig2, [virtual_right(2, :)], Rs(2) * scale_factor, 0.5, '-');
PlotCircle(fig2, [virtual_right(3, :)], Rs(3) * scale_factor, 0.5, '-');

vec = [cos(th + 1.75 * dtheta), sin(th + 1.75 * dtheta)];
facc = 0.97;
plabels = [(facc * R0) * vec; facc * (R0 + 1 * dx) * vec; facc * (R0 + 2 * dx) * vec; facc * (R0 + 3 * dx) * vec];
% Add lebels to the volumes
text(virtual_central(1, 1), virtual_central(1, 2), '$V_0$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);
text(virtual_central(2, 1), virtual_central(2, 2), '$V_1$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);
text(virtual_central(3, 1), virtual_central(3, 2), '$V_2$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);

% Add the labels along the curves
text(plabels(1, 1), plabels(1, 2), '$R_0$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);
text(plabels(2, 1), plabels(2, 2), '$R_1$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);
text(plabels(3, 1), plabels(3, 2), '$R_2$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);
text(plabels(4, 1), plabels(4, 2), '$R_3$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);

% Add circular sector for dxwall
PlotCircularSector(fig2, [0, 0], R0, th - dtheta / 2, th + dtheta / 2, 4, '-', 'r');
th_offset = +0.25 * dtheta;
dxlabelpos = 0.97 * R0 * [cos(th + th_offset), sin(th + th_offset)];
text(dxlabelpos(1), dxlabelpos(2), '$\Delta x$', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 14);
% plot marker particles
plot(pcentral(1), pcentral(2), 'ob', 'MarkerFaceColor', 'b');
plot(pright(1), pright(2), 'ob', 'MarkerFaceColor', 'b');
plot(pleft(1), pleft(2), 'ob', 'MarkerFaceColor', 'b');
axis([-2 * dx, 2 * dx, -R0 - 4 * dx, -R0 + 1.0 * dx]);
axis equal;
% set(gca, 'Visible', 'off');
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
% print -dpdf -painters epsFig
print(gcf, 'concave.pdf', '-dpdf', '-bestfit');

%%%% kernel function
h = 1;
rfine = linspace(0, 3 * h, 1000);
W = zeros(length(rfine), 1);
gradW = zeros(length(rfine), 1);

for k = 1:length(rfine)
    W(k) = kernel(rfine(k), h);
    gradW(k) = kernelGradient(rfine(k), h);
end

figure;
plot(rfine, W, 'k-'); hold on;
plot(rfine, gradW, 'r-');
xlabel('$r/h$');
legend('$W(r)$', '$\nabla W(r)$');
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
% print -dpdf -painters epsFig
print(gcf, 'kernel.pdf', '-dpdf', '-bestfit');

function W = kernel(r, h)

    q = r / h;

    K = 7.0/478.0 / pi / h / h; % 2d normalization

    tmp3 = (3.0 - q) .* (3.0 - q) .* (3.0 - q) .* (3.0 - q) .* (3.0 - q);
    tmp2 = (2.0 - q) .* (2.0 - q) .* (2.0 - q) .* (2.0 - q) .* (2.0 - q);
    tmp1 = (1.0 - q) .* (1.0 - q) .* (1.0 - q) .* (1.0 - q) .* (1.0 - q);

    if (q >= 0.0 & q < 1.0)
        W = K * (tmp3 - 6.0 * tmp2 + 15.0 * tmp1);
    elseif (q >= 1.0 & q < 2.0)
        W = K * (tmp3 - 6.0 * tmp2);
    elseif (q >= 2.0 & q < 3.0)
        W = K * tmp3;
    else
        W = 0.0;
    end

end

function W = kernelGradient(r, h)

    q = r / h;

    K = -5.0 * 7.0/478.0 / pi / h / h / h; % 2d normalization

    tmp3 = (3.0 - q) .* (3.0 - q) .* (3.0 - q) .* (3.0 - q);
    tmp2 = (2.0 - q) .* (2.0 - q) .* (2.0 - q) .* (2.0 - q);
    tmp1 = (1.0 - q) .* (1.0 - q) .* (1.0 - q) .* (1.0 - q);

    if (q >= 0.0 & q < 1.0)
        W = K * (tmp3 - 6.0 * tmp2 + 15.0 * tmp1);
    elseif (q >= 1.0 & q < 2.0)
        W = K * (tmp3 - 6.0 * tmp2);
    elseif (q >= 2.0 & q < 3.0)
        W = K * tmp3;
    else
        W = 0.0;
    end

end

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

function PlotCircle(fig, centre, radius, thickness, style)
    Nsamples = 1000;
    theta = linspace(0, 2 * pi, Nsamples);
    x = centre(1) + radius * cos(theta);
    y = centre(2) + radius * sin(theta);
    figure(fig);
    plot(x, y, 'k', 'LineWidth', thickness, 'LineStyle', style);
end

function PlotCircularSector(fig, centre, radius, theta0, thetaN, thickness, style, color)
    Nsamples = 1000;
    theta = linspace(theta0, thetaN, Nsamples);
    x = centre(1) + radius * cos(theta);
    y = centre(2) + radius * sin(theta);
    figure(fig);
    plot(x, y, 'LineWidth', thickness, 'LineStyle', style, 'Color', color);
end
