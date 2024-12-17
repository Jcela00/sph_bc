clear all; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 2);
set(groot, 'defaultFigureUnits', 'centimeters');
set(groot, 'defaultFigurePosition', [100, 100, 16.0, 10.0]); %double column

Rflower = 1;
Repicycloid = Rflower * 3;
Rhypocycloid = Repicycloid / 2;

kepi = 3;

[xepi, yepi] = Epicycloid(Repicycloid, kepi);
[xhyp, yhyp] = Hypocycloid(Rhypocycloid, kepi);

angles = [0:1:2] * 2 * pi / kepi + pi / kepi;
radial_dist = (5/3) * Repicycloid / 2;

figure; hold on; axis equal;
plot(xepi, yepi, 'r');
plot(xhyp, yhyp, 'b');

for n = 1:3
    centre = [radial_dist * cos(angles(n)), radial_dist * sin(angles(n))];
    [xstar, ystar] = PolarStarBurst(centre, Rflower, 0.4, 5, 7);
    plot(xstar, ystar, 'k');
end

axis tight;
set(gca, 'Visible', 'off');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, 'LatexFigures/EpicycloidGeometry.pdf', 'ContentType', 'vector', 'Resolution', 300);

% % plot radial sectors at angles

% for i = 1:length(angles)
%     plot([0, Repicycloid * cos(angles(i))], [0, Repicycloid * sin(angles(i))], 'r--');
%     plot([0, Rhypocycloid * cos(angles(i))], [0, Rhypocycloid * sin(angles(i))], 'b--');
% end

function [x y] = Epicycloid(R, k)
    t = linspace(0, 2 * pi, 1000);
    r = R / k;
    x = r * (k + 1) * cos(t) - r * cos((k + 1) * t);
    y = r * (k + 1) * sin(t) - r * sin((k + 1) * t);
end

function [x y] = Hypocycloid(R, k)
    t = linspace(0, 2 * pi, 1000);
    r = R / k;
    x = r * (k - 1) * cos(t) + r * cos((k - 1) * t);
    y = r * (k - 1) * sin(t) - r * sin((k - 1) * t);
end

function [x y] = PolarStarBurst(centre, a, b, k, m)
    t = linspace(0, 2 * pi, 1000);
    r = a * (1 + b * cos(k * t) .^ m);
    x = centre(1) + r .* cos(t);
    y = centre(2) + r .* sin(t);
end
