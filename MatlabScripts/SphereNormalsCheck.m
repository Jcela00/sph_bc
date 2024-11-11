clear all; clc; close all;

data = ParticleData(['Sphere_new_summation_ComputedNormals_30_30_80_2rf_1prc/file'], 0, ['sphere'], 1, 3, 1);

pos = data.Position(data.ObstacleIndexes, :);
normals = data.Normal(data.ObstacleIndexes, :);

[nobstacle, dim] = size(pos);
exact_normals = zeros(nobstacle, dim);
centre = data.ForceTransport(data.ObstacleIndexes, :);
centre = centre(1, :);

error = zeros(nobstacle, 1);

for k = 1:nobstacle
    exact_normals(k, :) = pos(k, :) - centre;
    exact_normals(k, :) = exact_normals(k, :) / norm(exact_normals(k, :));

    error(k) = norm(exact_normals(k, :) - normals(k, :));
end

figure; hold on;
scatter3(pos(:, 1), pos(:, 2), pos(:, 3), 20, 'filled', 'k', 'Displayname', 'Particles');
quiver3(pos(:, 1), pos(:, 2), pos(:, 3), normals(:, 1), normals(:, 2), normals(:, 3), 1, 'LineWidth', 2, 'Color', 'k', 'Displayname', 'Computed normals');
quiver3(pos(:, 1), pos(:, 2), pos(:, 3), exact_normals(:, 1), exact_normals(:, 2), exact_normals(:, 3), 1, 'LineWidth', 2, 'Color', 'r', 'Displayname', 'Exact normals');
axis equal;
legend;

average_error = mean(error)
max_error = max(error)
figure;
plot(error * 100, 'LineWidth', 2);
ylabel('Error %');
xlabel('particle index');
yline(average_error * 100, 'k', 'Average error', 'LineWidth', 2);
