clear all; clc; close all;
% Load the data in csv format
data = csvread('channel.csv', 1, 0);
% delete nan density values
data = data(~isnan(data(:, 2)), :);
type_particle = data(:, 1);
density = data(:, 2);
pressure = data(:, 3);
velocity = data(:, 4:6);
domain = data(:, 7);
positions = data(:, 8:10);

Nparticles = length(type_particle);
% sph parameters
dp = 0.0085;
H = sqrt(3) * dp;
mass = 0.000614125;

plot3(positions(:, 1), positions(:, 2), positions(:, 3), 'o');
% compute limits in coordinates
max_x = max(positions(:, 1));
min_x = min(positions(:, 1));
max_y = max(positions(:, 2));
min_y = min(positions(:, 2));
max_z = max(positions(:, 3));
min_z = min(positions(:, 3));

% positons for profile computation
x_samples = 80;
x_grid = linspace(min_x, max_x, x_samples);
y_grid = min_y * ones(1, x_samples); % y is constant (2D)

z_samples = 1;
z_points = linspace(min_z, max_z, z_samples + 4);
z_points = z_points(3:end - 2);

velocity_interp = zeros(x_samples, 3, z_samples);

for j = 1:z_samples % loop over z coordinates along the channel
    zval = z_points(j);
    z_grid = zval * ones(1, x_samples);
    R_grid = [x_grid; y_grid; z_grid]';

    for i = 1:x_samples % loop over all x coordinates across the channel

        for k = 1:Nparticles % loop over all particles

            if type_particle(k) == 1 % fluid particles
                r = norm(R_grid(i, :) - positions(k, :));
                W = kernel(r, H);

                velocity_interp(i, :, j) = velocity_interp(i, :, j) + mass * velocity(k, :) * W / density(k, 1);

            end

        end

    end

end

x_fine = linspace(min_x + 2.5 * dp, max_x - 2.5 * dp, 1000);

for j = 1:z_samples
    umax = max(-velocity_interp(:, 3, j));
    w_a = prof_a(x_fine, umax, [min_x, max_x]);
    figure;
    plot(x_grid, -velocity_interp(:, 3, j), 'k'); hold on;
    plot(x_fine, w_a, 'r--')
    xlabel('x'); ylabel('Velocity z');
    title(['z = ', num2str(z_points(j))]);
    xline(min_x + 2.5 * dp, 'b--');
    xline(max_x - 2.5 * dp, 'b--');
    particle_initial = min_x:dp:max_x;
    plot(particle_initial, zeros(1, length(particle_initial)), 'bo');
end

% figure;
% plot(x_grid, velocity_grid(:, 2), 'k'); xlabel('x'); ylabel('Velocity y');

% figure;
% plot(x_grid, velocity_grid(:, 3), 'k'); xlabel('x'); ylabel('Velocity z');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytical velocity profile

function u = prof_a(x, umax, xlims)
    % u =g/(2*nu)*x*(h-x)
    % or umax= g*h^2/(8*nu) and then
    % u = (4umax/h)*x(1-x/h)

    xmin = xlims(1);
    xmax = xlims(2);

    dp = 0.0085;
    x_in_min = xmin + 2.5 * dp;
    x_in_max = xmax - 2.5 * dp;

    g = 9.81;
    % nu = 0.001;
    h = 0.17;
    % w = (g / (2 * nu)) .* x .* (h - x);
    u = zeros(1, length(x));

    u = (4 * umax / h) * x .* (1 - x / h);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = kernel(r, h)
    % Kernel function
    q = r / h;
    K = 1 / (pi * h ^ 3);
    %K = 10 / (7 * pi * h ^ 2);

    if q < 1
        W = K * (1 - 1.5 * q ^ 2 + 0.75 * q ^ 3);
    elseif q < 2
        W = K * 0.25 * (2 - q) ^ 3;
    else
        W = 0;
    end

end
