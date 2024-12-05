classdef ParticleData

    properties
        % CSV PROPERTIES TO BE LOADED
        Type
        Density
        Pressure
        DeltaDensity
        Velocity
        VelocityTransport
        Force
        ForceTransport
        Normal
        Volumes
        Omega
        Vorticity
        Domain
        Position
        Time
        % Curvature % Not used anymore
        % ArcLength
        %%%%%%%%%%%%
        filename_base
        BaseName
        PlotName
        dim
        VelMax
        VelMin
        Npart
        NpartFluid
        mass
        dp
        H
        rho0
        Hconst_

        BoundaryIndexes
        ObstacleIndexes
        FluidIndexes

        TypeBoundary = 0;
        TypeObstacle = 1;
        TypeFluid = 2;
    end

    methods

        function obj = ParticleData(filename, filenum, titlename, rho0, dim, Hconst)
            % Load the data in csv format
            obj.filename_base = filename;
            name = [filename '_' num2str(filenum) '.csv'];
            data = csvread(name, 1, 0);

            obj.BaseName = titlename;
            obj.PlotName = [obj.BaseName]; % ' ' num2str(filenum)];
            obj.dim = dim;
            obj.Hconst_ = Hconst;
            obj.rho0 = rho0;
            obj.Time = data(1, 1);
            data = data(:, 2:end); % remove time column
            obj.Type = data(:, 1);
            obj.Density = data(:, 2);
            obj.Pressure = data(:, 3);
            obj.DeltaDensity = data(:, 4);
            obj.Velocity = data(:, 5:7);
            obj.VelocityTransport = data(:, 8:10);
            obj.Force = data(:, 11:13);
            obj.ForceTransport = data(:, 14:16);
            obj.Normal = data(:, 17:19);
            obj.Volumes = data(:, 20:22);
            obj.Omega = data(:, 23);
            obj.Vorticity = data(:, 24);
            % skip columns red_vel red_fx red_fy 25 26 27
            obj.Domain = data(:, 28);
            obj.Position = data(:, 29:31);

            obj.Npart = length(obj.Type);
            obj.NpartFluid = sum(obj.Type == obj.TypeFluid);
            obj.VelMax = max(sqrt(sum(obj.Velocity .^ 2, 2)));
            obj.VelMin = min(sqrt(sum(obj.Velocity .^ 2, 2)));

            %Indexes
            obj.BoundaryIndexes = find(obj.Type == obj.TypeBoundary);
            obj.ObstacleIndexes = find(obj.Type == obj.TypeObstacle);
            obj.FluidIndexes = find(obj.Type == obj.TypeFluid);

            obj.dp = obj.Volumes(obj.FluidIndexes(1), 1); % fluid particles contain the dp value in vol[1]

            if (dim == 2)
                obj.mass = rho0 * obj.dp ^ 2;
            elseif (dim == 3)
                obj.mass = rho0 * obj.dp ^ 3;
            else
                error('Invalid dimension parameter, must be 2 or 3.');
            end

            obj.H = Hconst * obj.dp;

        end

        function obj = Update(obj, filenum)
            % create a new instance of the particle data object at a new time step
            obj = ParticleData(obj.filename_base, filenum, obj.BaseName, obj.rho0, obj.dim, obj.Hconst_);
        end

        function obj = OnlyBoundaryFluid(obj, fluidorboundary)

            if (fluidorboundary ~= 0 && fluidorboundary ~= 1)
                error('Invalid input, must be 0 or 1.');
            end

            typetmp = obj.Type;

            obj.Type = typetmp(typetmp == fluidorboundary);
            obj.Density = obj.Density(typetmp == fluidorboundary);
            obj.Pressure = obj.Pressure(typetmp == fluidorboundary);
            obj.DeltaDensity = obj.DeltaDensity(typetmp == fluidorboundary);
            obj.Velocity = obj.Velocity(typetmp == fluidorboundary, :);
            obj.VelocityTransport = obj.VelocityTransport(typetmp == fluidorboundary, :);
            obj.Force = obj.Force(typetmp == fluidorboundary, :);
            obj.ForceTransport = obj.ForceTransport(typetmp == fluidorboundary, :);
            obj.Normal = obj.Normal(typetmp == fluidorboundary, :);
            obj.Volumes = obj.Volumes(typetmp == fluidorboundary, :);
            obj.Omega = obj.Omega(typetmp == fluidorboundary);
            obj.Vorticity = obj.Vorticity(typetmp == fluidorboundary);
            obj.Domain = obj.Domain(typetmp == fluidorboundary);
            obj.Position = obj.Position(typetmp == fluidorboundary, :);
            obj.Npart = length(obj.Type);

        end

        function [obj zmatrix_cell] = GridDataInterpVelocity(obj, xi, yi, H)

            x = obj.Position(obj.FluidIndexes, 1);
            y = obj.Position(obj.FluidIndexes, 2);
            u = obj.Velocity(obj.FluidIndexes, 1);
            v = obj.Velocity(obj.FluidIndexes, 2);

            Z = obj.Velocity(:, 1);

            [xq, yq] = meshgrid(xi, yi);

            ui = griddata(x, y, u, xq, yq, 'cubic');
            vi = griddata(x, y, v, xq, yq, 'cubic');

            zmatrix_cell = {ui, vi};
        end

        function [obj, zmatrix_cell] = SPHinterpVelocity(obj, xi, yi, H, setting)
            % magnitude is the property to be interpolated
            % function to interpolate particle properties to a grid defined by the cartesian product of xi and yi
            % the interpolation is done using the SPH kernel considering particles inside the kernel support
            % the result is returned in zmatrix

            % setting = 0 for x velocity interpolation
            % setting = 1 for y velocity interpolation
            % setting = 2 for both
            % setting = 3 for velocity magnitude

            if (setting == 0 || setting == 1)

                magnitude = obj.Velocity(:, setting + 1);
                zmatrix = zeros(length(xi), length(yi));

                total = length(xi) * length(yi);

                for i = 1:length(xi)

                    for j = 1:length(yi)

                        x = xi(i); y = yi(j);
                        Wsum = 0;

                        for k = 1:length(obj.FluidIndexes) % loop over fluid particles
                            n = obj.FluidIndexes(k); % acces the fluid particle index

                            r = norm([x, y] - obj.Position(n, 1:2));

                            if (r < 3 * H)
                                W = obj.kernel(r, H);
                                Wsum = Wsum + W * obj.mass / obj.Density(n, 1);
                                zmatrix(i, j) = zmatrix(i, j) + obj.mass * magnitude(n, 1) * W / obj.Density(n, 1);
                            end

                        end

                        zmatrix(i, j) = zmatrix(i, j) / Wsum;

                        progress = ((i - 1) * length(yi) + j) / total * 100;
                        fprintf('\rInterpolating: %.2f%%', progress);
                    end

                end

                zmatrix_cell = {zmatrix};

            elseif (setting == 2)
                zmatrix_u = zeros(length(xi), length(yi));
                zmatrix_v = zeros(length(xi), length(yi));
                zmatrix_w = zeros(length(xi), length(yi));

                pos = obj.Position(obj.FluidIndexes, 1:2);
                vel = obj.Velocity(obj.FluidIndexes, 1:2);
                vorticity = obj.Vorticity(obj.FluidIndexes, :);
                density = obj.Density(obj.FluidIndexes, 1);
                total = length(xi) * length(yi);

                for i = 1:length(xi)

                    for j = 1:length(yi)

                        x = xi(i); y = yi(j);
                        Wsum = 0;

                        for k = 1:length(obj.FluidIndexes) % loop over fluid particles
                            % n = obj.FluidIndexes(k); % acces the fluid particle index

                            r = norm([x, y] - pos(k, :));

                            if (r < 3 * H)
                                W = obj.kernel(r, H);
                                Wsum = Wsum + W * obj.mass / density(k, 1);
                                zmatrix_u(i, j) = zmatrix_u(i, j) + obj.mass * vel(k, 1) * W / density(k, 1);
                                zmatrix_v(i, j) = zmatrix_v(i, j) + obj.mass * vel(k, 2) * W / density(k, 1);
                                zmatrix_w(i, j) = zmatrix_w(i, j) + obj.mass * vorticity(k, 1) * W / density(k, 1);

                            end

                        end

                        zmatrix_u(i, j) = zmatrix_u(i, j) / Wsum;
                        zmatrix_v(i, j) = zmatrix_v(i, j) / Wsum;
                        zmatrix_w(i, j) = zmatrix_w(i, j) / Wsum;

                        progress = ((i - 1) * length(yi) + j) / total * 100;
                        fprintf('\rInterpolating: %.2f%%', progress);
                    end

                end

                zmatrix_cell = {zmatrix_u, zmatrix_v, zmatrix_w};

            elseif (setting == 3)
                magnitude = sqrt(obj.Velocity(:, 1) .^ 2 + obj.Velocity(:, 2) .^ 2);
                zmatrix = zeros(length(xi), length(yi));

                total = length(xi) * length(yi);

                for i = 1:length(xi)

                    for j = 1:length(yi)

                        x = xi(i); y = yi(j);
                        Wsum = 0;

                        for k = 1:length(obj.FluidIndexes) % loop over fluid particles
                            n = obj.FluidIndexes(k); % acces the fluid particle index

                            r = norm([x, y] - obj.Position(n, 1:2));

                            if (r < 3 * H)
                                W = obj.kernel(r, H);
                                Wsum = Wsum + W * obj.mass / obj.Density(n, 1);
                                zmatrix(i, j) = zmatrix(i, j) + obj.mass * magnitude(n, 1) * W / obj.Density(n, 1);
                            end

                        end

                        zmatrix(i, j) = zmatrix(i, j) / Wsum;

                        progress = ((i - 1) * length(yi) + j) / total * 100;
                        fprintf('\rInterpolating: %.2f%%', progress);
                    end

                end

                zmatrix_cell = {zmatrix};
            else
                error('Invalid setting parameter, must be 0, 1, 2 or 3.');
            end

        end

        function [obj, zmatrix_cell] = SPHinterpVelocityWithObstacle(obj, xi, yi, H)
            % magnitude is the property to be interpolated
            % function to interpolate particle properties to a grid defined by the cartesian product of xi and yi
            % the interpolation is done using the SPH kernel considering particles inside the kernel support
            % the result is returned in zmatrix

            zmatrix_u = zeros(length(xi), length(yi));
            zmatrix_v = zeros(length(xi), length(yi));
            zmatrix_w = zeros(length(xi), length(yi));

            pos = obj.Position(obj.FluidIndexes, 1:2);
            vel = obj.Velocity(obj.FluidIndexes, 1:2);
            vorticity = obj.Vorticity(obj.FluidIndexes, :);
            density = obj.Density(obj.FluidIndexes, 1);
            total = length(xi) * length(yi);

            % need to skip the interpolation points inside obstacle.
            obstacle_centre = obj.ForceTransport(obj.ObstacleIndexes, 1:2);
            obstacle_centre = obstacle_centre(1, :);
            square_half_side = 0.5;

            for i = 1:length(xi)

                for j = 1:length(yi)

                    x = xi(i); y = yi(j);

                    [~, isInside] = obj.IsInsideSquare([x, y], obstacle_centre, square_half_side);

                    if (isInside)
                        continue;
                    end

                    Wsum = 0;

                    for k = 1:length(obj.FluidIndexes) % loop over fluid particles
                        % n = obj.FluidIndexes(k); % acces the fluid particle index

                        r = norm([x, y] - pos(k, :));

                        if (r < 3 * H)
                            W = obj.kernel(r, H);
                            Wsum = Wsum + W * obj.mass / density(k, 1);
                            zmatrix_u(i, j) = zmatrix_u(i, j) + obj.mass * vel(k, 1) * W / density(k, 1);
                            zmatrix_v(i, j) = zmatrix_v(i, j) + obj.mass * vel(k, 2) * W / density(k, 1);
                            zmatrix_w(i, j) = zmatrix_w(i, j) + obj.mass * vorticity(k, 1) * W / density(k, 1);

                        end

                    end

                    zmatrix_u(i, j) = zmatrix_u(i, j) / Wsum;
                    zmatrix_v(i, j) = zmatrix_v(i, j) / Wsum;
                    zmatrix_w(i, j) = zmatrix_w(i, j) / Wsum;

                    progress = ((i - 1) * length(yi) + j) / total * 100;
                    fprintf('\rInterpolating: %.2f%%', progress);
                end

            end

            zmatrix_cell = {zmatrix_u, zmatrix_v, zmatrix_w};

        end

        function [obj, isInside] = IsInsideSquare(obj, point, centre, half_side)

            cx = centre(1); cy = centre(2);
            x = point(1); y = point(2);
            % Check if the point is within the square boundaries
            isInside = (x >= cx - half_side) && (x <= cx + half_side) && (y >= cy - half_side) && (y <= cy + half_side);
        end

        function [obj, isInside] = IsInsideTriangle(obj, point, centre, sideLength)

            % Extract point and center coordinates
            px = point(1);
            py = point(2);
            centre_x = centre(1);
            centre_y = centre(2);

            left_vertex = [centre_x - sqrt(3) * sideLength / 3, centre_y];
            top_vertex = [centre_x + sqrt(3) * sideLength / 6, centre_y + sideLength / 2];
            bottom_vertex = [centre_x + sqrt(3) * sideLength / 6, centre_y - sideLength / 2];

            % Extract vertex coordinates
            x1 = left_vertex(1); y1 = left_vertex(2);
            x2 = top_vertex(1); y2 = top_vertex(2);
            x3 = bottom_vertex(1); y3 = bottom_vertex(2);

            % Compute barycentric coordinates
            denominator = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
            a = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / denominator;
            b = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / denominator;
            c = 1 - a - b;

            % Check if the point is inside the triangle
            isInside = (a >= 0) && (b >= 0) && (c >= 0);

        end

        % function [obj, zmatrix_cell] = SPHinterpVorticity(obj, xi, yi, H)
        %     % function to interpolate particle properties to a grid defined by the cartesian product of xi and yi
        %     % the interpolation is done using the SPH kernel considering particles inside the kernel support
        %     % the result is returned in zmatrix

        %     pos = obj.Position(obj.FluidIndexes, 1:2);
        %     vorticity = obj.Vorticity(obj.FluidIndexes, :);
        %     density = obj.Density(obj.FluidIndexes, 1);
        %     mass = obj.mass;

        %     zmatrix = zeros(length(xi), length(yi));

        %     total = length(xi) * length(yi);

        %     for i = 1:length(xi)

        %         for j = 1:length(yi)

        %             x = xi(i); y = yi(j);
        %             Wsum = 0;

        %             for k = 1:length(obj.FluidIndexes) % loop over fluid particles

        %                 r = norm([x, y] - pos(k, :));

        %                 if (r < 3 * H)
        %                     W = obj.kernel(r, H);
        %                     Wsum = Wsum + W * mass / density(k, 1);
        %                     zmatrix(i, j) = zmatrix(i, j) + mass * vorticity(k, 1) * W / density(k, 1);
        %                 end

        %             end

        %             zmatrix(i, j) = zmatrix(i, j) / Wsum;

        %             progress = ((i - 1) * length(yi) + j) / total * 100;
        %             fprintf('\rInterpolating Vorticity: %.2f%%', progress);
        %         end

        %     end

        %     zmatrix_cell = {zmatrix};

        % end

        function obj = PlotVelocityContour(obj, fig, levels, lnwdth)
            X = obj.Position(:, 1);
            Y = obj.Position(:, 2);
            Z = obj.Velocity(:, 1);

            N = 50;
            xi = linspace(0, 0.1, N);
            yi = linspace(0, 0.1, N);

            Vmag = sqrt(obj.Velocity(:, 1) .^ 2 + obj.Velocity(:, 2) .^ 2);
            [obj, Zi] = obj.SPHinterpVelocity(xi, yi, obj.H, 3);

            [xi, yi] = meshgrid(xi, yi);
            Zi = Zi';
            figure(fig);
            contour(xi, yi, Zi, levels, "LineWidth", lnwdth, "LineColor", "black", 'ShowText', 'on');

        end

        function obj = PlotStreamlinesOnly(obj, fig2, lnwdth, xrange, yrange, hfac, Npoints, levelsStreamfunction)

            N = sqrt(length(obj.FluidIndexes));

            Nx = Npoints(1);
            Ny = Npoints(2);

            Hfac = hfac;
            % Hfac = N / Npoints;

            % equispace interpolation points
            xi = linspace(xrange(1), xrange(2), Nx);
            yi = linspace(yrange(1), yrange(2), Ny);

            fprintf('Interpolating...\n');
            [obj, v_cell] = obj.SPHinterpVelocity(xi, yi, Hfac * obj.H, 2);

            ui = v_cell{1};
            vi = v_cell{2};
            vorticity = v_cell{3};

            % add boundary velues to ui and vi
            % u(x=0) = 0    u(x=1) = 0    u(y=0) = 0    u(y=1) = 1
            % v(x=0) = 0    v(x=1) = 0    v(y=0) = 0    v(y=1) = 0

            vorticity = vorticity';

            dx = xi(2) - xi(1);
            dy = yi(2) - yi(1);

            [xi, yi] = meshgrid(xi, yi);

            size(vorticity)

            % SOLVE STREAM FUNCTION
            psi = zeros(Ny, Nx); % Initial guess for stream function
            tol = 1e-10; % Convergence tolerance
            maxIter = 20000; % Maximum number of iterations
            err = 1; iter = 0;
            % remove NaN values in vorticity
            % vorticity(isnan(vorticity)) = 0;

            while err > tol && iter < maxIter
                psi_old = psi;

                for i = 2:Ny - 1

                    for j = 2:Nx - 1

                        % Solve using Jacobi method
                        if (~isnan(vorticity(i, j)))

                            psi(i, j) = 0.25 * (psi_old(i + 1, j) + psi_old(i - 1, j) + psi_old(i, j + 1) + psi_old(i, j - 1) + dx ^ 2 * vorticity(i, j));

                        end

                    end

                end

                % Compute error as the maximum change
                err = max(max(abs(psi - psi_old)));
                iter = iter + 1;
            end

            fprintf('Jacobi solver in %d iterations with error %e\n', iter, err);

            % % normalize the vectors for quiver plot
            % norm = sqrt(ui .^ 2 + vi .^ 2);
            % unorm = ui ./ norm;
            % vnorm = vi ./ norm;

            % save('streamlines.mat', 'xi', 'yi', 'ui', 'vi', 'vorticity');

            % figure(fig1);
            % l = streamslice(xi, yi, ui, vi, streamlinedensity, 'noarrows', 'cubic');
            % set(l, 'Color', 'k', 'LineWidth', lnwdth);
            % axis equal;
            % axis([0 1 0 1]);
            % xline(0, 'k'); yline(0, 'k');
            % xline(1, 'k'); yline(1, 'k');
            % xlabel('$x$'); ylabel('$y$');
            % set(gca, 'FontSize', 11);
            % set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);

            figure(fig2);
            contour(xi, yi, psi, 50, 'LineWidth', lnwdth, 'LineColor', 'black');
            axis equal;
            axis([xrange(1) xrange(2) yrange(1) yrange(2)]);
            xlabel('$x$'); ylabel('$y$');
            set(gca, 'FontSize', 11);
            set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);

            % figure(fig2);
            % contour(xi, yi, -vorticity, levelsVorticity, 'ShowText', 'on', 'LineWidth', lnwdth, 'LineColor', 'black', 'LabelSpacing', 10000);
            % axis equal;
            % axis([0 1 0 1]);
            % xlabel('$x$'); ylabel('$y$');
            % set(gca, 'FontSize', 11);
            % set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);

            % figure(fig3);
            % % Cavity.PlotParticles(fig3, mrksz);
            % contourf(xi, yi, sqrt(ui .^ 2 + vi .^ 2), 100, 'LineStyle', 'none');
            % % quiver(xi, yi, unorm, vnorm, 'Color', 'm', 'LineWidth', 1);
            % axis equal;
            % fig3.Colormap = turbo;
            % cb = colorbar;
            % set(cb, 'TickLabelInterpreter', 'latex');
            % clim([obj.VelMin, obj.VelMax]);
            % axis([0 1 0 1]);
            % xlabel('$x$'); ylabel('$y$');
            % set(gca, 'FontSize', 11);
            % set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);

        end

        function obj = PlotStreamlines(obj, fig2, fig3, fig4, lnwdth, xrange, yrange, Npoints, streamlinedensity, levelsStreamfunction, levelsVorticity)
            % this function was intended for the cavity and it got too complex, there is a separate plotstreamlinesonly just for streamlines
            N = sqrt(length(obj.FluidIndexes));

            Hfac = N / Npoints

            % equispace interpolation points
            xi = linspace(xrange(1), xrange(2), Npoints);
            yi = linspace(yrange(1), yrange(2), Npoints);

            fprintf('Interpolating...\n');
            [obj, v_cell] = obj.SPHinterpVelocity(xi, yi, Hfac * obj.H, 2);
            % [obj, v_cell] = obj.GridDataInterpVelocity(xi, yi, obj.H);

            ui = v_cell{1};
            vi = v_cell{2};
            vorticity = v_cell{3};

            % add boundary velues to ui and vi
            % u(x=0) = 0    u(x=1) = 0    u(y=0) = 0    u(y=1) = 1
            % v(x=0) = 0    v(x=1) = 0    v(y=0) = 0    v(y=1) = 0

            % [nrows, ncols] = size(vi);
            % vi = [zeros(1, ncols); vi; zeros(1, ncols)];
            % vi = [zeros(nrows + 2, 1) vi zeros(nrows + 2, 1)];
            % ui = [zeros(1, ncols); ui; zeros(1, ncols)];
            % ui = [zeros(nrows + 2, 1) ui ones(nrows + 2, 1)];

            vi(1, :) = zeros(1, Npoints);
            vi(end, :) = zeros(1, Npoints);
            vi(:, 1) = zeros(Npoints, 1);
            vi(:, end) = zeros(Npoints, 1);

            ui(1, :) = zeros(1, Npoints);
            ui(end, :) = zeros(1, Npoints);
            ui(:, 1) = zeros(Npoints, 1);
            ui(:, end) = ones(Npoints, 1);

            ui = ui';
            vi = vi';

            vorticity = vorticity';

            dx = xi(2) - xi(1);
            dy = yi(2) - yi(1);

            [xi, yi] = meshgrid(xi, yi);

            % SOLVE STREAM FUNCTION
            psi = zeros(Npoints, Npoints); % Initial guess for stream function
            tol = 1e-10; % Convergence tolerance
            maxIter = 20000; % Maximum number of iterations
            err = 1; iter = 0;

            while err > tol && iter < maxIter
                psi_old = psi;

                for i = 2:Npoints - 1

                    for j = 2:Npoints - 1
                        % Solve using Jacobi method
                        psi(i, j) = 0.25 * (psi_old(i + 1, j) + psi_old(i - 1, j) + psi_old(i, j + 1) + psi_old(i, j - 1) + dx ^ 2 * vorticity(i, j));
                    end

                end

                % Compute error as the maximum change
                err = max(max(abs(psi - psi_old)));
                iter = iter + 1;
            end

            fprintf('Jacobi solver in %d iterations with error %e\n', iter, err);

            % normalize the vectors for quiver plot
            norm = sqrt(ui .^ 2 + vi .^ 2);
            unorm = ui ./ norm;
            vnorm = vi ./ norm;

            % save('streamlines.mat', 'xi', 'yi', 'ui', 'vi', 'vorticity');

            % figure(fig1);
            % l = streamslice(xi, yi, ui, vi, streamlinedensity, 'noarrows', 'cubic');
            % set(l, 'Color', 'k', 'LineWidth', lnwdth);
            % axis equal;
            % axis([0 1 0 1]);
            % xline(0, 'k'); yline(0, 'k');
            % xline(1, 'k'); yline(1, 'k');
            % xlabel('$x$'); ylabel('$y$');
            % set(gca, 'FontSize', 11);
            % set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);

            figure(fig4);
            contour(xi, yi, psi, levelsStreamfunction, 'LineWidth', lnwdth, 'LineColor', 'black');
            axis equal;
            axis([0 1 0 1]);
            xlabel('$x$'); ylabel('$y$');
            set(gca, 'FontSize', 11);
            set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);

            figure(fig2);
            contour(xi, yi, -vorticity, levelsVorticity, 'ShowText', 'on', 'LineWidth', lnwdth, 'LineColor', 'black', 'LabelSpacing', 10000);
            axis equal;
            axis([0 1 0 1]);
            xlabel('$x$'); ylabel('$y$');
            set(gca, 'FontSize', 11);
            set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);

            figure(fig3);
            % Cavity.PlotParticles(fig3, mrksz);
            contourf(xi, yi, sqrt(ui .^ 2 + vi .^ 2), 100, 'LineStyle', 'none');
            % quiver(xi, yi, unorm, vnorm, 'Color', 'm', 'LineWidth', 1);
            axis equal;
            fig3.Colormap = turbo;
            cb = colorbar;
            set(cb, 'TickLabelInterpreter', 'latex');
            clim([obj.VelMin, obj.VelMax]);
            axis([0 1 0 1]);
            xlabel('$x$'); ylabel('$y$');
            set(gca, 'FontSize', 11);
            set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);

        end

        function plotNormals(obj, fig)

            figure(fig); hold on;

            for k = 1:length(obj.Position(:, 1))

                if (obj.Type(k) == obj.TypeBoundary || obj.Type(k) == obj.TypeObstacle)
                    plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'y', 'Marker', 'o', 'MarkerSize', 6);
                    % quiver(obj.Position(k, 1), obj.Position(k, 2), obj.Normal(k, 1), obj.Normal(k, 2), obj.Nfluid(1) * obj.dp * 0.025, 'Color', 'k');
                    % plot([obj.Position(k, 1) obj.Position(k, 1) - 2.5 * obj.dp * obj.Normal(k, 1)], [obj.Position(k, 2) obj.Position(k, 2) - 2.5 * obj.dp * obj.Normal(k, 2)], '--k');
                    % plot(obj.Position(k, 1) -0.5 * obj.dp * obj.Normal(k, 1), obj.Position(k, 2) -0.5 * obj.dp * obj.Normal(k, 2), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 's', 'MarkerSize', 6);
                    % plot(obj.Position(k, 1) -1.5 * obj.dp * obj.Normal(k, 1), obj.Position(k, 2) -1.5 * obj.dp * obj.Normal(k, 2), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 's', 'MarkerSize', 6);
                    % plot(obj.Position(k, 1) -2.5 * obj.dp * obj.Normal(k, 1), obj.Position(k, 2) -2.5 * obj.dp * obj.Normal(k, 2), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 's', 'MarkerSize', 6);
                end

            end

        end

        function plotNormals3d(obj, fig)

            figure(fig); hold on;

            for k = 1:length(obj.Position(:, 1))

                if (obj.Type(k) == obj.TypeBoundary || obj.Type(k) == obj.TypeObstacle)
                    plot(obj.Position(k, 1), obj.Position(k, 2), obj.Position(k, 3), 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'y', 'Marker', 'o', 'MarkerSize', 6);
                    quiver(obj.Position(k, 1), obj.Position(k, 2), obj.Position(k, 3), obj.Normal(k, 1), obj.Normal(k, 2), obj.Normal(k, 3), obj.Nfluid(1) * obj.dp * 0.025, 'Color', 'k');
                    % plot([obj.Position(k, 1) obj.Position(k, 1) - 2.5 * obj.dp * obj.Normal(k, 1)], [obj.Position(k, 2) obj.Position(k, 2) - 2.5 * obj.dp * obj.Normal(k, 2)], '--k');
                    % plot(obj.Position(k, 1) -0.5 * obj.dp * obj.Normal(k, 1), obj.Position(k, 2) -0.5 * obj.dp * obj.Normal(k, 2), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 's', 'MarkerSize', 6);
                    % plot(obj.Position(k, 1) -1.5 * obj.dp * obj.Normal(k, 1), obj.Position(k, 2) -1.5 * obj.dp * obj.Normal(k, 2), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 's', 'MarkerSize', 6);
                    % plot(obj.Position(k, 1) -2.5 * obj.dp * obj.Normal(k, 1), obj.Position(k, 2) -2.5 * obj.dp * obj.Normal(k, 2), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 's', 'MarkerSize', 6);
                end

            end

        end

        function PlotParticles(obj, fig, mrksz)
            figure(fig);
            custom_colormap = turbo;

            for k = 1:obj.Npart

                if (obj.Type(k) == obj.TypeFluid)
                    vel_magnitude = norm(obj.Velocity(k, :));
                    color_index = round(1 + ((255 - 1) / (obj.VelMax - obj.VelMin)) * (vel_magnitude - obj.VelMin));
                    plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', custom_colormap(color_index, :), 'MarkerFaceColor', custom_colormap(color_index, :), 'Marker', 'o', 'MarkerSize', mrksz);
                else
                    plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Marker', 'o', 'MarkerSize', mrksz);

                end

            end

            fig.Colormap = custom_colormap;

            cb = colorbar;
            set(cb, 'TickLabelInterpreter', 'latex');
            clim([obj.VelMin, obj.VelMax]);
            % title(obj.PlotName);

        end

        function PlotParticlesVorticity(obj, fig, mrksz, cmap, wlimits)
            figure(fig);
            custom_colormap = cmap;
            wmax = wlimits(2);
            wmin = wlimits(1);

            for k = 1:obj.Npart

                if (obj.Type(k) == obj.TypeFluid)
                    magnitude = obj.Vorticity(k);
                    color_index = round(1 + ((size(custom_colormap, 1) - 1) / (wmax - wmin)) * (magnitude - wmin));
                    color_index = max(1, min(color_index, size(custom_colormap, 1))); % Clamp to [1, 256]
                    plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', custom_colormap(color_index, :), 'MarkerFaceColor', custom_colormap(color_index, :), 'Marker', 'o', 'MarkerSize', mrksz);
                else
                    plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Marker', 'o', 'MarkerSize', mrksz);

                end

                progress = k / obj.Npart * 100;
                fprintf('\rPlotting Particles: %.2f%%', progress);

            end

            fig.Colormap = custom_colormap;

            cb = colorbar;
            set(cb, 'TickLabelInterpreter', 'latex');
            clim([wmin, wmax]);
            % title(obj.PlotName);

        end

        function PlotParticlesRange(obj, fig, mrksz, xrange, yrange)
            figure(fig);

            % m = 256;
            % control_points = [0.0, 0.231, 0.298, 0.752;
            %                   0.5, 1.000, 1.000, 1.000;
            %                   1.0, 0.706, 0.016, 0.150; ];

            % x = control_points(:, 1);
            % r = control_points(:, 2);
            % g = control_points(:, 3);
            % b = control_points(:, 4);
            % cmap = [interp1(x, r, linspace(0, 1, m))', interp1(x, g, linspace(0, 1, m))', interp1(x, b, linspace(0, 1, m))'];

            cmap = obj.BlueGreenYellowRed(256);
            custom_colormap = cmap;

            indexInside = obj.IsInRange(xrange, yrange);

            obj.Npart = length(indexInside);
            obj.Type = obj.Type(indexInside);
            obj.Position = obj.Position(indexInside, :);
            obj.Velocity = obj.Velocity(indexInside, :);
            obj.VelMax = max(sqrt(obj.Velocity(:, 1) .^ 2 + obj.Velocity(:, 2) .^ 2));
            obj.VelMin = 0;

            for k = 1:obj.Npart

                if (obj.Type(k) == obj.TypeFluid)
                    vel_magnitude = norm(obj.Velocity(k, :));
                    color_index = round(1 + ((255 - 1) / (obj.VelMax - obj.VelMin)) * (vel_magnitude - obj.VelMin));
                    plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', custom_colormap(color_index, :), 'MarkerFaceColor', custom_colormap(color_index, :), 'Marker', 'o', 'MarkerSize', mrksz);
                else
                    plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Marker', 'o', 'MarkerSize', mrksz);

                end

            end

            fig.Colormap = custom_colormap;

            cb = colorbar;
            set(cb, 'TickLabelInterpreter', 'latex');
            clim([obj.VelMin, obj.VelMax]);
            % title(obj.PlotName);

        end

        function PlotParticlesVorticityRange(obj, fig, mrksz, wrange, xrange, yrange)
            figure(fig);

            cmap = obj.paraview_red_blue_colormap(256);
            custom_colormap = cmap;
            wmin = wrange(1);
            wmax = wrange(2);

            indexInside = obj.IsInRange(xrange, yrange);

            obj.Npart = length(indexInside);
            obj.Type = obj.Type(indexInside);
            obj.Position = obj.Position(indexInside, :);
            obj.Vorticity = obj.Vorticity(indexInside, :);

            for k = 1:obj.Npart

                if (obj.Type(k) == obj.TypeFluid)
                    magnitude = obj.Vorticity(k);
                    color_index = round(1 + ((size(custom_colormap, 1) - 1) / (wmax - wmin)) * (magnitude - wmin));
                    color_index = max(1, min(color_index, size(custom_colormap, 1))); % Clamp to [1, 256]
                    plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', custom_colormap(color_index, :), 'MarkerFaceColor', custom_colormap(color_index, :), 'Marker', 'o', 'MarkerSize', mrksz);
                else
                    plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Marker', 'o', 'MarkerSize', mrksz);

                end

                progress = k / obj.Npart * 100;
                fprintf('\rPlotting Particles: %.2f%%', progress);

            end

            fig.Colormap = custom_colormap;

            cb = colorbar('southoutside');
            set(cb, 'TickLabelInterpreter', 'latex');
            clim([wmin, wmax]);
            % title(obj.PlotName);

        end

        function cmap = paraview_red_blue_colormap(obj, n)

            white_range = [-0.4, 0.4]; % Default white range

            % Define the RGB values for the minimum, central, and maximum colors
            rgb_min = [59, 76, 192] / 255; % Minimum blue
            rgb_mid = [222, 220, 218] / 255; % Central ghost white
            rgb_max = [180, 4, 38] / 255; % Maximum red

            % Define normalized positions for the color interpolation
            x = linspace(-1, 1, n); % Normalized colormap range [-1, 1]
            white_min = white_range(1); % Start of white range
            white_max = white_range(2); % End of white range

            % Initialize the colormap
            cmap = zeros(n, 3);

            % Create the colormap by interpolating RGB values
            for i = 1:3
                % Blue to white
                idx_blue = x <= white_min;
                cmap(idx_blue, i) = interp1([-1, white_min], [rgb_min(i), rgb_mid(i)], x(idx_blue));

                % White region
                idx_white = (x > white_min) & (x < white_max);
                cmap(idx_white, i) = rgb_mid(i);

                % White to red
                idx_red = x >= white_max;
                cmap(idx_red, i) = interp1([white_max, 1], [rgb_mid(i), rgb_max(i)], x(idx_red));
            end

        end

        function cmap = BlueGreenYellowRed(obj, nColors)
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

        function cmap = BlueYellowRedCmap(obj, nColors)
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

        function [indexInside] = IsInRange(obj, xrange, yrange)
            indexInside = [];

            for k = 1:obj.Npart

                if (obj.Position(k, 1) >= xrange(1) && obj.Position(k, 1) <= xrange(2) && obj.Position(k, 2) >= yrange(1) && obj.Position(k, 2) <= yrange(2))
                    indexInside = [indexInside; k];
                end

            end

        end

        function PlotCurvatureNormals(obj, fig, mrksz, indexsetting, CurvatureLimits, vectorsize)

            m = 256;

            % Define the RGB control points for ParaView's Red-Blue colormap
            control_points = [
                              0.0, 0.231, 0.298, 0.752; % Dark blue
                              0.5, 1.000, 1.000, 1.000; % White
                              1.0, 0.706, 0.016, 0.150; % Dark red
                              ];

            % Interpolate the colormap
            x = control_points(:, 1); % Position of control points
            r = control_points(:, 2); % Red values
            g = control_points(:, 3); % Green values
            b = control_points(:, 4); % Blue values

            % Interpolate to create the colormap with 'm' colors
            cmap = [interp1(x, r, linspace(0, 1, m))', interp1(x, g, linspace(0, 1, m))', interp1(x, b, linspace(0, 1, m))'];

            figure(fig);
            custom_colormap = cmap;

            if (indexsetting == 0)
                indexes = obj.ObstacleIndexes;
            elseif (indexsetting == 1)
                indexes = obj.BoundaryIndexes;
            end

            position = obj.Position(indexes, 1:2);
            curvature = obj.Volumes(indexes, 2);
            normal = obj.Normal(indexes, 1:2);

            if obj.BaseName == "ConcaveCylinder"
                % need to filter out particles with r < 0.75
                position = obj.Position(indexes, 1:2);
                r = sqrt(position(:, 1) .^ 2 + position(:, 2) .^ 2);
                indexes = indexes(r > 0.75);
                position = obj.Position(indexes, 1:2);
                curvature = obj.Volumes(indexes, 2);
                normal = obj.Normal(indexes, 1:2);
            end

            dx = obj.Volumes(1, 1);
            fprintf('name: %s, dx: %f\n', obj.BaseName, dx);
            maxcurv = max(curvature);
            mincurv = min(curvature);
            fprintf('Measured maxcurv: %f, mincurv: %f\n', maxcurv, mincurv);

            % if (maxcurv == mincurv)
            %     maxcurv = maxcurv + 1;
            %     mincurv = mincurv - 1;
            % end

            % obj.BaseName
            % mincurv = min(curvature)
            % maxcurv = max(curvature)

            mincurv = CurvatureLimits(1);
            maxcurv = CurvatureLimits(2);

            for k = 1:length(indexes)

                curv = curvature(k);
                pos = position(k, :);
                norm = normal(k, :);

                color_index = round(1 + ((255 - 1) / (maxcurv - mincurv)) * (curv - mincurv));
                plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', custom_colormap(color_index, :), 'Marker', 'o', 'MarkerSize', mrksz);
                quiver(pos(1), pos(2), norm(1), norm(2), vectorsize, 'Color', 'k', 'LineWidth', 2);

            end

            fig.Colormap = custom_colormap;
            cb = colorbar;
            set(cb, 'TickLabelInterpreter', 'latex');
            clim([mincurv, maxcurv]);
            set(gca, 'Visible', 'off');
            set(gca, 'FontSize', 11);
            set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
            % title(obj.PlotName);

        end

        function W = kernel(obj, r, H)

            if nargin < 3
                H = obj.H;
                disp('Using default H');
            end

            q = r / H;

            if (obj.dim == 2)
                K = 7.0/478.0 / pi / H / H;
            elseif (obj.dim == 3)
                K = 1.0/120.0 / pi / H / H / H;
            end

            tmp3 = (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q);
            tmp2 = (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q);
            tmp1 = (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q);

            if (q >= 0.0 && q < 1.0)
                W = K * (tmp3 - 6.0 * tmp2 + 15.0 * tmp1);
            elseif (q >= 1.0 && q < 2.0)
                W = K * (tmp3 - 6.0 * tmp2);
            elseif (q >= 2.0 && q < 3.0)
                W = K * tmp3;
            else
                W = 0.0;
            end

        end

        function CheckSum = CheckKernel(obj)

            CheckSum = zeros(obj.Npart, 1);

            figure; hold on;

            for k = 1:obj.Npart

                % only consider fluid particles
                if (obj.Type(k, 1) == obj.TypeFluid)

                    % Ignore particles near the periodic boundary
                    if ((obj.Position(k, 3) > obj.PositionLimits(3, 1) + 2 * obj.H) && (obj.Position(k, 3) < obj.PositionLimits(3, 2) - 2 * obj.H))

                        ra = obj.Position(k, :);

                        for j = 1:obj.Npart
                            rb = obj.Position(j, :);
                            r = norm(ra - rb);

                            if (r < 2 * obj.H)
                                W = obj.kernel(r);
                                CheckSum(k, 1) = CheckSum(k, 1) + W * obj.mass / obj.Density(j, 1);

                            end

                        end

                        plot(k, CheckSum(k, 1), 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
                    end

                end

            end

            % remove zero values
            CheckSum = CheckSum(CheckSum ~= 0);
        end

        % function plotMagnitude(obj, magnitude, name, fig, kk)

        %     figure(fig); hold on;

        %     if length(magnitude(1, :)) == 1 % scalar
        %         % do nothing
        %     else
        %         magnitude = sum(magnitude, 2);
        %     end

        %     for k = 1:obj.Npart

        %         if magnitude(k) == 0
        %             plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'Marker', 'o', 'MarkerSize', 6);
        %         else
        %             plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'o', 'MarkerSize', 10);
        %         end

        %     end

        %     xlabel('$x$'); ylabel('$y$');
        %     title([name num2str(kk)]);

        % end

    end

end
