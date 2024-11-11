clear all; clc; close all;
load('streamlines.mat');

% test to compute the streamfunction from the interpolated vorticity field
% there are two options, take the interpolated vorticity as rhs or compute the vorticity
% from the interpolated velocity, both seem to produce the same result

[nx, ny] = size(vorticity);
x = linspace(0, 1, nx);
y = linspace(0, 1, ny);
dx = x(2) - x(1);
dy = y(2) - y(1);

omega = vorticity;

% [X, Y] = meshgrid(x, y);
% % u, v: 2D arrays of velocity components (given data)
% omega = zeros(nx, ny);
% u = ui; v = vi;

% [dvdx, ~] = gradient(v, dx, dy);
% [~, dudy] = gradient(u, dx, dy);
% omega = dvdx - dudy;

psi = zeros(nx, ny); % Initial guess for stream function
tol = 1e-10; % Convergence tolerance
maxIter = 20000; % Maximum number of iterations
err = 1; iter = 0;

while err > tol && iter < maxIter
    psi_old = psi;

    for i = 2:nx - 1

        for j = 2:ny - 1
            psi(i, j) = 0.25 * (psi_old(i + 1, j) + psi_old(i - 1, j) + ...
                psi_old(i, j + 1) + psi_old(i, j - 1) + ...
                dx ^ 2 * omega(i, j));
        end

    end

    % Compute error as the maximum change
    err = max(max(abs(psi - psi_old)));
    iter = iter + 1;
end

fprintf('Converged in %d iterations with error %e\n', iter, err);

levels = [-1e-10 -1e-7 -1e-5 -1e-4 -0.01 -0.03 -0.05 -0.07 -0.09 -0.1 -0.11 -0.115 -0.1175 1e-8 1e-7 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];
contour(x, y, psi, levels, 'LineWidth', 1.0, 'Color', 'k');
axis equal;
