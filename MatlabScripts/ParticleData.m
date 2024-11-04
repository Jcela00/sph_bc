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

        function [obj, zmatrix] = SPHinterp(obj, xi, yi)
            %function to interpolate particle properties to a grid defined by the cartesian product of xi and yi
            % the interpolation is done using the SPH kernel considering particles inside the kernel support
            % the result is returned in zmatrix

            zmatrix = zeros(length(xi), length(yi));

            for i = 1:length(xi)

                for j = 1:length(yi)
                    x = xi(i);
                    y = yi(j);
                    Wsum = 0;

                    for k = 1:length(obj.FluidIndexes) % loop over fluid particles
                        % acces the fluid particle index
                        n = FluidIndexes(k);

                        if (obj.Type(n) == obj.TypeFluid)
                            r = norm([x, y] - obj.Position(n, 1:2));

                            if (r < 3 * obj.H)
                                W = obj.kernel(r);
                                Wsum = Wsum + W * obj.mass / obj.Density(n, 1);
                                zmatrix(i, j) = zmatrix(i, j) + obj.mass * obj.Velocity(n, 1) * W / obj.Density(n, 1);
                            end

                        end

                    end

                    zmatrix(i, j) = zmatrix(i, j) / Wsum;

                end

            end

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

        function obj = PlotContour(obj, fig, levels)
            X = obj.Position(:, 1);
            Y = obj.Position(:, 2);
            Z = obj.Velocity(:, 1);

            N = round(1.5 * sqrt(obj.NpartFluid));
            xi = linspace(0, 0.1, N);
            yi = linspace(0, 0.1, N);
            % [xi, yi] = meshgrid(xi, yi);
            % Zi = griddata(X, Y, Z, xi, yi, 'cubic');

            [obj, Zi] = obj.SPHinterp(xi, yi);

            [xi, yi] = meshgrid(xi, yi);
            Zi = Zi';
            figure(fig);
            contour(xi, yi, Zi, levels, "LineWidth", 4, "LineColor", "black", 'ShowText', 'on');

            xlabel('$x$');
            ylabel('$y$');
            title(obj.PlotName);
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

        function PlotParticles(obj, fig)
            figure(fig);
            custom_colormap = turbo;

            for k = 1:obj.Npart

                if (obj.Type(k) == obj.TypeFluid)
                    vel_magnitude = norm(obj.Velocity(k, :));
                    color_index = round(1 + ((255 - 1) / (obj.VelMax - obj.VelMin)) * (vel_magnitude - obj.VelMin));
                    plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', custom_colormap(color_index, :), 'MarkerFaceColor', custom_colormap(color_index, :), 'Marker', 'o', 'MarkerSize', 4);
                else
                    plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'Marker', 'o', 'MarkerSize', 4);

                end

            end

            fig.Colormap = custom_colormap;

            colorbar;
            clim([obj.VelMin, obj.VelMax]);
            xlabel('$z$'); ylabel('$x$');
            title(obj.PlotName);

        end

        function W = kernel(obj, r)

            q = r / obj.H;

            if (obj.dim == 2)
                K = 7.0/478.0 / pi / obj.H / obj.H;
            elseif (obj.dim == 3)
                K = 1.0/120.0 / pi / obj.H / obj.H / obj.H;
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
