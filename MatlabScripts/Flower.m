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

% a = 1;
% b = 0.3;
% k = 5;
% m = 5;

% [x, y] = RoseCurve(a, b, k, m);

% figure; hold on;
% plot(x, y, 'k', 'DisplayName', 'Rose Curve');
% axis equal;
% title('Rose Curve');

% t = linspace(0, 2 * pi, 100000);
% dt = t(2) - t(1);

% accumulated_arc = 0;

% for n = 1:length(t)
%     R = RoseCurveRadius(a, b, k, m, t(n));
%     DR = RoseCurveRadiusDerivative(a, b, k, m, t(n));
%     ds = sqrt(R ^ 2 + DR ^ 2) * dt;
%     accumulated_arc = accumulated_arc + ds;
% end

% fprintf('Total arc length: %f\n', accumulated_arc);
% Npoints = 120;
% arc_element = accumulated_arc / Npoints;

% t = linspace(0, 2 * pi, 1e7);
% dt = t(2) - t(1);

% accumulated_arc = 0;
% particles = zeros(2, Npoints);
% % add first point
% R = RoseCurveRadius(a, b, k, m, t(1));
% particles(:, 1) = [R * cos(t(1)); R * sin(t(1))];
% count = 2;

% RPrev = R;
% DrPrev = RoseCurveRadiusDerivative(a, b, k, m, t(1));

% f_prev = sqrt(RPrev ^ 2 + DrPrev ^ 2);

% for n = 1:length(t)
%     R = RoseCurveRadius(a, b, k, m, t(n));
%     DR = RoseCurveRadiusDerivative(a, b, k, m, t(n));
%     f = sqrt(R ^ 2 + DR ^ 2);

%     ds = (f + f_prev) * dt / 2;
%     f_prev = f;

%     accumulated_arc = accumulated_arc + ds;

%     if accumulated_arc >= arc_element
%         accumulated_arc = 0;
%         particles(:, count) = [R * cos(t(n)); R * sin(t(n))];
%         count = count + 1;
%     end

%     % if count > Npoints
%     %     break;
%     % end

% end

% % % remove zero columns from particles
% % particles = particles(:, 1:count - 1);

% plot(particles(1, :), particles(2, :), 'ro', 'DisplayName', 'Particles', 'MarkerSize', 10);
% legend('Location', 'best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

figure; hold on;
k = 3;
R = 1;
Nsamples = 1e7;

t = linspace(0, 2 * pi, Nsamples);
dt = t(2) - t(1);

accumulated_arc_epi = 0;
accumulated_arc_hipo = 0;

for n = 1:length(t)
    [dx1 dy1] = EpicycloidDerivative(1, k, t(n));
    [dx2 dy2] = HypocycloidDerivative(0.5, k, t(n));

    ds1 = sqrt(dx1 ^ 2 + dy1 ^ 2) * dt;
    ds2 = sqrt(dx2 ^ 2 + dy2 ^ 2) * dt;

    accumulated_arc_epi = accumulated_arc_epi + ds1;
    accumulated_arc_hipo = accumulated_arc_hipo + ds2;
end

fprintf('Epicycloid arc length: %f\n', accumulated_arc_epi);
fprintf('Hypocycloid arc length: %f\n', accumulated_arc_hipo);

% discretize Epicycloid

Npoints = 120;
arc_element = accumulated_arc_epi / Npoints;

t = linspace(0, 2 * pi, Nsamples);
dt = t(2) - t(1);

accumulated_arc = 0;
particles = zeros(2, Npoints);
% add first point
[x y] = Epicycloid(R, k, 0);

particles(:, 1) = [x; y];

count = 2;

[dx dy] = EpicycloidDerivative(R, k, 0);
f_prev = sqrt(dx ^ 2 + dy ^ 2);

for n = 1:length(t)
    [dx dy] = EpicycloidDerivative(R, k, t(n));

    f = sqrt(dx ^ 2 + dy ^ 2);

    ds = (f + f_prev) * dt / 2;
    f_prev = f;

    accumulated_arc = accumulated_arc + ds;

    if accumulated_arc >= arc_element
        accumulated_arc = 0;
        [x y] = Epicycloid(R, k, t);
        particles(:, count) = [x(n); y(n)];
        count = count + 1;
    end

end

% % remove zero columns from particles
% particles = particles(:, 1:count - 1);

plot(particles(1, :), particles(2, :), 'ro', 'DisplayName', 'Particles', 'MarkerSize', 10);

% discretize Hypocycloid

R = 0.5;
k = 3;

Npoints = 120;
arc_element = accumulated_arc_hipo / Npoints;

t = linspace(0, 2 * pi, Nsamples);
dt = t(2) - t(1);

accumulated_arc = 0;
particles = zeros(2, Npoints);

% add first point
[x y] = Hypocycloid(R, k, 0);

particles(:, 1) = [x; y];

count = 2;

[dx dy] = HypocycloidDerivative(R, k, 0);
f_prev = sqrt(dx ^ 2 + dy ^ 2);

for n = 1:length(t)
    [dx dy] = HypocycloidDerivative(R, k, t(n));

    f = sqrt(dx ^ 2 + dy ^ 2);

    ds = (f + f_prev) * dt / 2;
    f_prev = f;

    accumulated_arc = accumulated_arc + ds;

    if accumulated_arc >= arc_element
        accumulated_arc = 0;
        [x y] = Hypocycloid(R, k, t);
        particles(:, count) = [x(n); y(n)];
        count = count + 1;
    end

end

plot(particles(1, :), particles(2, :), 'ro', 'DisplayName', 'Particles', 'MarkerSize', 10);
legend('Location', 'best');
axis equal;

% close all;

% figure; hold on;
% t = linspace(0, 2 * pi, 100);

% for n = 1:length(t)
%     [x y] = Epicycloid(R, k, t(n));
%     plot(x, y, 'or');
%     [x y] = Hypocycloid(R, k, t(n));
%     plot(x, y, 'ob');
%     title(['t = ', num2str(t(n))]);
%     pause(0.1)
% end

function [x y] = RoseCurve(a, b, k, m)
    t = linspace(0, 2 * pi, 1000);
    r = a * (1 + b * cos(k * t) .^ m);
    x = r .* cos(t);
    y = r .* sin(t);
end

function R = RoseCurveRadius(a, b, k, m, t)
    R = a * (1 + b * cos(k * t) .^ m);
end

function DR = RoseCurveRadiusDerivative(a, b, k, m, t)
    DR = -a * b * m * k * cos(k * t) .^ (m - 1) .* sin(k * t);
end

function [x y] = Epicycloid(R, k, t)
    r = R / k;
    x = r * (k + 1) * cos(t) - r * cos((k + 1) * t);
    y = r * (k + 1) * sin(t) - r * sin((k + 1) * t);
end

function [dx dy] = EpicycloidDerivative(R, k, t)
    r = R / k;
    dx = -r * (k + 1) * sin(t) + r * (k + 1) * sin((k + 1) * t);
    dy = r * (k + 1) * cos(t) - r * (k + 1) * cos((k + 1) * t);
end

function [x y] = Hypocycloid(R, k, t)
    r = R / k;
    x = r * (k - 1) * cos(t) + r * cos((k - 1) * t);
    y = r * (k - 1) * sin(t) - r * sin((k - 1) * t);
end

function [dx dy] = HypocycloidDerivative(R, k, t)
    r = R / k;
    dx = -r * (k - 1) * sin(t) - r * (k - 1) * sin((k - 1) * t);
    dy = r * (k - 1) * cos(t) - r * (k - 1) * cos((k - 1) * t);
end
