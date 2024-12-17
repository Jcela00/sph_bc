clear all; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 2);
% This is for exportgraphics
set(groot, 'defaultFigureUnits', 'centimeters');
% set(groot, 'defaultFigurePosition', [0, 0, 8.5, 6.0]); %single column
set(groot, 'defaultFigurePosition', [100, 100, 16.0, 10.0]); %double column

Repi = 1;
k = 3;
[x, y] = Epicycloid(Repi, k);
[nx ny] = EpicycloidNormal(Repi, k);

figure; hold on;
plot(x, y, 'r');
% quiver(x, y, nx, ny, 'k');
axis equal;
Rhipo = 1;
[x, y] = Hypocycloid(Rhipo, k);
[nx ny] = HypocycloidNormals(Rhipo, k);

plot(x, y, 'b');
% quiver(x, y, nx, ny, 'k');
axis equal;

n = 0:k - 1
maxangles = pi / k + 2 * pi * n / k
[xmax ymax] = Epicycloid_t(Repi, 3, maxangles(1));
dmax = (sqrt(xmax .^ 2 + ymax .^ 2));
dmax = dmax

% plot radial lines at these angles
for i = 1:length(maxangles)
    x = [0 dmax * cos(maxangles(i))];
    y = [0 dmax * sin(maxangles(i))];
    plot(x, y, 'k');
end

title(['Epicycloid and Hypocycloid k=' num2str(k)]);

figure; hold on;
k = 11;
epsilon = 0.1;
[x, y] = PolygonalApproximation(epsilon, k);
plot(x, y, 'k');
axis equal;
title(['Polygonal Approximation k=' num2str(k)]);

figure; hold on;
epsilon = 0.5;
[x, y] = Treefoil(epsilon);
plot(x, y, 'k', 'DisplayName', 'Treefoil');
axis equal;
title('Treefoil');

figure; hold on;
[x, y] = Astroid();
plot(x, y, 'k', 'DisplayName', 'Astroid');
axis equal;
title('Astroid');

figure; hold on;
a = 1;
k = 6;
epsilon = 0.7;
[x, y] = RoundedStar(a, k, epsilon);
plot(x, y, 'k');
axis equal;
title(['Rounded Star k=' num2str(k)]);

figure; hold on;
k = 3;
n = 1.01;
[x, y] = StarflowerSuperellipse(k, n);
plot(x, y, 'k', 'DisplayName', ['n = ', num2str(n)]);
n = 6;
[x, y] = StarflowerSuperellipse(k, n);
plot(x, y, 'b', 'DisplayName', ['n = ', num2str(n)]);
legend('Location', 'best');
axis equal;
title('Starflower Superellipse');

figure; hold on;
a = 1;
b = 0.5;
k = 5;
m = 101;

[x, y] = PolarStarBurst(a, b, k, m);
plot(x, y, 'k');
axis equal;
title(['Generalized Polar Star k = ' num2str(k) ' m = ' num2str(m)]);

figure; hold on;
epsilon = 0.2;
p = 20;
[x, y] = HighFreqPerturbedCircle(epsilon, p);
plot(x, y, 'k', 'DisplayName', 'High Frequency Perturbed Circle');
axis equal;
title('High Frequency Perturbed Circle');

figure; hold on;
[x, y] = FlowerOfBohr();
plot(x, y, 'k', 'DisplayName', 'Flower of Bohr');
axis equal;
title('Flower of Bohr');

function [x y] = Epicycloid(R, k)
    t = linspace(0, 2 * pi, 1000);
    r = R / k;
    x = r * (k + 1) * cos(t) - r * cos((k + 1) * t);
    y = r * (k + 1) * sin(t) - r * sin((k + 1) * t);
end

function [nx ny] = EpicycloidNormal(R, k)
    t = linspace(0, 2 * pi, 1000);
    r = R / k;
    dx = -r * (k + 1) * sin(t) + r * (k + 1) * sin((k + 1) * t);
    dy = r * (k + 1) * cos(t) - r * (k + 1) * cos((k + 1) * t);
    nx = -dy;
    ny = dx;
end

function [x y] = Epicycloid_t(R, k, t)
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

function [nx ny] = HypocycloidNormals(R, k)
    t = linspace(0, 2 * pi, 1000);
    r = R / k;
    dx = -r * (k - 1) * sin(t) - r * (k - 1) * sin((k - 1) * t);
    dy = r * (k - 1) * cos(t) - r * (k - 1) * cos((k - 1) * t);
    nx = dy;
    ny = -dx;
end

function [x y] = Hypotrochoid(R, r, d)
    t = linspace(0, 10 * pi, 1000);
    x = (R - r) * cos(t) + d * cos((R - r) / r * t);
    y = (R - r) * sin(t) - d * sin((R - r) / r * t);
end

function [x y] = PolygonalApproximation(epsilon, k)
    t = linspace(0, 2 * pi, 1000);
    R = 1 + epsilon * cos(k * t);
    x = R .* cos(t);
    y = R .* sin(t);
end

function [x y] = Treefoil(epsilon)
    t = linspace(0, 2 * pi, 1000);
    r = 1 + epsilon * cos(3 * t);
    x = r .* cos(t);
    y = r .* sin(t);
end

function [x y] = Astroid()
    t = linspace(0, 2 * pi, 1000);
    x = cos(t) .^ 3;
    y = sin(t) .^ 3;
end

function [x y] = RoundedStar(a, k, epsilon)
    t = linspace(0, 2 * pi, 1000);
    r = a * (1 + epsilon * cos(k * t));
    x = r .* cos(t);
    y = r .* sin(t);
end

function [x y] = StarflowerSuperellipse(k, n)
    t = linspace(0, 2 * pi, 1000);
    r = (abs(cos(k * t)) .^ n + abs(sin(k * t)) .^ n) .^ (-1 / n);
    x = r .* cos(t);
    y = r .* sin(t);
end

function [x y] = PolarStarBurst(a, b, k, m)
    t = linspace(0, 2 * pi, 1000);
    r = a * (1 + b * cos(k * t) .^ m);
    x = r .* cos(t);
    y = r .* sin(t);
end

function [x y] = HighFreqPerturbedCircle(epsilon, p)
    t = linspace(0, 2 * pi, 1000);
    r = 1 + epsilon * cos(p * t);
    x = r .* cos(t);
    y = r .* sin(t);
end

function [x y] = FlowerOfBohr()
    t = linspace(0, 2 * pi, 1000);
    x = cos(t) + 0.35 * cos(6 * t);
    y = sin(t) + 0.35 * sin(6 * t);
end
