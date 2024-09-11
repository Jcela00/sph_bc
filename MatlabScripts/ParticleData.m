classdef ParticleData

    properties
        BaseName
        PlotName
        dim
        Type
        Density
        Pressure
        DeltaDensity
        Force
        Velocity
        VelMax
        VelMin
        ForceTransport
        VelocityTransport
        Normal
        Curvature
        ArcLength
        Domain
        Volumes
        Omega
        Position
        Npart
        NpartFluid
        mass
        dp
        H
        rho0
        Nfluid
        Nboundary
        PositionLimits
        FluidBox
        filename_base
        dim_self
        Hconst_self
    end

    methods

        function obj = ParticleData(filename, filenum, titlename, dp, Nfluid, Nboundary, rho0, dim, Hconst)
            % Load the data in csv format
            obj.filename_base = filename;
            name = [filename '_' num2str(filenum) '.csv'];
            data = csvread(name, 1, 0);

            % delete nan density values
            % data = data(~isnan(data(:, 2)), :);
            obj.BaseName = titlename;
            obj.PlotName = [obj.BaseName]; % ' ' num2str(filenum)];
            obj.dim = dim;
            obj.Type = data(:, 1);
            obj.Density = data(:, 2);
            obj.Pressure = data(:, 3);
            obj.DeltaDensity = data(:, 4);
            obj.Force = data(:, 5:7);
            obj.Velocity = data(:, 8:10);
            obj.VelMax = max(sqrt(sum(obj.Velocity .^ 2, 2)));
            obj.VelMin = min(sqrt(sum(obj.Velocity .^ 2, 2)));
            obj.ForceTransport = data(:, 11:13);
            obj.VelocityTransport = data(:, 14:16);
            obj.Normal = data(:, 17:19);
            obj.Curvature = data(:, 20);
            obj.ArcLength = data(:, 21);
            obj.Volumes = data(:, 22:24);
            obj.Omega = data(:, 25);
            obj.Domain = data(:, 26);
            obj.Position = data(:, 27:29);
            obj.Npart = length(obj.Type);
            obj.dp = dp;
            obj.dim_self = dim;
            obj.Hconst_self = Hconst;
            obj.NpartFluid = sum(obj.Type == 1);

            if (dim == 2)
                obj.mass = rho0 * obj.dp ^ 2;
            elseif (dim == 3)
                obj.mass = rho0 * obj.dp ^ 3;
            else
                error('Invalid dimension parameter, must be 2 or 3.');
            end

            obj.H = Hconst * obj.dp;
            obj.rho0 = rho0;
            obj.Nfluid = Nfluid;
            obj.Nboundary = Nboundary;
            obj.PositionLimits = [min(obj.Position(:, 1)) max(obj.Position(:, 1));
                                  min(obj.Position(:, 2)) max(obj.Position(:, 2));
                                  min(obj.Position(:, 3)) max(obj.Position(:, 3))];
            obj.FluidBox = [Nfluid(1) * dp
                            Nfluid(2) * dp
                            Nfluid(3) * dp];

        end

        function [obj, zmatrix] = SPHinterp(obj, xi, yi)

            zmatrix = zeros(length(xi), length(yi));

            for i = 1:length(xi)

                for j = 1:length(yi)
                    x = xi(i);
                    y = yi(j);
                    Wsum = 0;

                    for n = 1:obj.Npart

                        if (obj.Type(n) == 1)
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
            obj.Force = obj.Force(typetmp == fluidorboundary, :);
            obj.Velocity = obj.Velocity(typetmp == fluidorboundary, :);
            obj.ForceTransport = obj.ForceTransport(typetmp == fluidorboundary, :);
            obj.VelocityTransport = obj.VelocityTransport(typetmp == fluidorboundary, :);
            obj.Normal = obj.Normal(typetmp == fluidorboundary, :);
            obj.Curvature = obj.Curvature(typetmp == fluidorboundary);
            obj.ArcLength = obj.ArcLength(typetmp == fluidorboundary);
            obj.Volumes = obj.Volumes(typetmp == fluidorboundary, :);
            obj.Omega = obj.Omega(typetmp == fluidorboundary);
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

        % function obj = ReorderData (obj, Isort)
        %     obj.Type = obj.Type(Isort);
        %     obj.Density = obj.Density(Isort);
        %     obj.Pressure = obj.Pressure(Isort);
        %     obj.DeltaDensity = obj.DeltaDensity(Isort);
        %     obj.Force = obj.Force(Isort, :);
        %     obj.Velocity = obj.Velocity(Isort, :);
        %     obj.ForceTransport = obj.ForceTransport(Isort, :);
        %     obj.VelocityTransport = obj.VelocityTransport(Isort, :);
        %     obj.Normal = obj.Normal(Isort, :);
        %     obj.Curvature = obj.Curvature(Isort);
        %     obj.ArcLength = obj.ArcLength(Isort);
        %     % obj.InteractCount = obj.InteractCount(Isort);
        %     obj.Domain = obj.Domain(Isort);
        %     obj.Position = obj.Position(Isort, :);
        % end

        % function [objDiff, pack] = CompDiff(obj1, obj2) % must be well sorted
        %     objDiff = obj1;
        %     objDiff.Type = abs(obj1.Type - obj2.Type);
        %     objDiff.Density = abs(obj1.Density - obj2.Density);
        %     objDiff.Pressure = abs(obj1.Pressure - obj2.Pressure);
        %     objDiff.DeltaDensity = abs(obj1.DeltaDensity - obj2.DeltaDensity);
        %     objDiff.Force = abs(obj1.Force - obj2.Force);
        %     objDiff.Velocity = abs(obj1.Velocity - obj2.Velocity);
        %     objDiff.ForceTransport = abs(obj1.ForceTransport - obj2.ForceTransport);
        %     objDiff.VelocityTransport = abs(obj1.VelocityTransport - obj2.VelocityTransport);
        %     objDiff.Normal = abs(obj1.Normal - obj2.Normal);
        %     % objDiff.InteractCount = abs(obj1.InteractCount - obj2.InteractCount);
        %     objDiff.Curvature = abs(obj1.Curvature - obj2.Curvature);
        %     objDiff.ArcLength = abs(obj1.ArcLength - obj2.ArcLength);

        %     totalE_rho = sum(objDiff.Density);
        %     totalE_p = sum(objDiff.Pressure);
        %     totalE_Force = sum(norm(objDiff.Force));
        %     totalE_Velocity = sum(norm(objDiff.Velocity));
        %     totalE_ForceTransport = sum(norm(objDiff.ForceTransport));
        %     totalE_VelocityTransport = sum(norm(objDiff.VelocityTransport));
        %     totalE_Normal = sum(norm(objDiff.Normal));
        %     totalE_Curvature = sum(objDiff.Curvature);
        %     totalE_ArcLength = sum(objDiff.ArcLength);
        %     totalE_InteractCount = sum(objDiff.InteractCount);
        %     number_error = abs(obj1.Npart - obj2.Npart);

        %     pack = [totalE_rho, totalE_p, totalE_Force, totalE_Velocity, totalE_ForceTransport, totalE_VelocityTransport, totalE_Normal, totalE_Curvature, totalE_ArcLength, totalE_InteractCount, number_error];

        % end

        function plotNormals(obj, fig)

            figure(fig); hold on;

            for k = 1:length(obj.Position(:, 1))

                if (obj.Type(k) == 0)
                    plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'y', 'Marker', 'o', 'MarkerSize', 6);
                    % quiver(obj.Position(k, 1), obj.Position(k, 2), obj.Normal(k, 1), obj.Normal(k, 2), obj.Nfluid(1) * obj.dp * 0.025, 'Color', 'k');
                    % plot([obj.Position(k, 1) obj.Position(k, 1) - 2.5 * obj.dp * obj.Normal(k, 1)], [obj.Position(k, 2) obj.Position(k, 2) - 2.5 * obj.dp * obj.Normal(k, 2)], '--k');
                    % plot(obj.Position(k, 1) -0.5 * obj.dp * obj.Normal(k, 1), obj.Position(k, 2) -0.5 * obj.dp * obj.Normal(k, 2), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 's', 'MarkerSize', 6);
                    % plot(obj.Position(k, 1) -1.5 * obj.dp * obj.Normal(k, 1), obj.Position(k, 2) -1.5 * obj.dp * obj.Normal(k, 2), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 's', 'MarkerSize', 6);
                    % plot(obj.Position(k, 1) -2.5 * obj.dp * obj.Normal(k, 1), obj.Position(k, 2) -2.5 * obj.dp * obj.Normal(k, 2), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 's', 'MarkerSize', 6);
                end

            end

        end

        function plotMagnitude(obj, magnitude, name, fig, kk)

            figure(fig); hold on;

            if length(magnitude(1, :)) == 1 % scalar
                % do nothing
            else
                magnitude = sum(magnitude, 2);
            end

            for k = 1:obj.Npart

                if magnitude(k) == 0
                    plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'Marker', 'o', 'MarkerSize', 6);
                else
                    plot(obj.Position(k, 1), obj.Position(k, 2), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'o', 'MarkerSize', 10);
                end

            end

            xlabel('$x$'); ylabel('$y$');
            title([name num2str(kk)]);

        end

        function obj = Update(obj, filenum)
            % create a new instance of the particle data object at a new time step
            obj = ParticleData(obj.filename_base, filenum, obj.BaseName, obj.dp, obj.Nfluid, obj.Nboundary, obj.rho0, obj.dim, obj.Hconst_self);
        end

        function PlotParticles(obj, fig)
            figure(fig);
            custom_colormap = turbo;

            for k = 1:obj.Npart

                if (obj.Type(k) == 1)
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

        function ScatterParticlePlot(obj, fig, colormarker, markersize)
            figure(fig);

            for k = 1:obj.Npart

                if obj.Type(k) == 1

                    if (k == 1)
                        plot(obj.Position(:, 1), obj.Velocity(:, 3), colormarker, 'MarkerSize', markersize, 'DisplayName', obj.PlotName);
                    else
                        plot(obj.Position(:, 1), obj.Velocity(:, 3), colormarker, 'MarkerSize', markersize, 'HandleVisibility', 'off');
                    end

                    % else

                    %     if (k == 1)
                    %         plot(obj.Position(:, 1), obj.VelocityTransport(:, 3), marker, 'MarkerSize', markersize, 'DisplayName', obj.PlotName);

                    %     else
                    %         plot(obj.Position(:, 1), obj.VelocityTransport(:, 3), marker, 'MarkerSize', markersize, 'HandleVisibility', 'off');
                    %     end

                end

            end

            xlabel('$x$'); ylabel('Velocity $z$');
            xline(obj.PositionLimits(1, 1), 'k', 'HandleVisibility', 'off');
            xline(obj.PositionLimits(1, 2), 'k--', 'HandleVisibility', 'off');
            xline(0.0, 'k', 'HandleVisibility', 'off');
            xline(obj.FluidBox(1), 'k', 'HandleVisibility', 'off');
            particle_initial = obj.PositionLimits(1, 1):obj.dp:obj.PositionLimits(1, 2);
            plot(particle_initial, zeros(1, length(particle_initial)), 'bo', 'HandleVisibility', 'off');

        end

        function W = kernelCubic(obj, r)
            % Kernel function
            q = r / obj.H;

            if (obj.dim == 2)
                K = 10 / (7 * pi * obj.H ^ 2);
            elseif (obj.dim == 3)
                K = 1 / (pi * obj.H ^ 3);
            end

            if q < 1
                W = K * (1 - 1.5 * q ^ 2 + 0.75 * q ^ 3);
            elseif q < 2
                W = K * 0.25 * (2 - q) ^ 3;
            else
                W = 0;
            end

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

        function PlotProfile(obj, fig, x_samples, z_coord)

            % Sample positons for profile computation
            x_sample_grid = linspace(obj.PositionLimits(1, 1), obj.PositionLimits(1, 2), x_samples);

            y_grid = obj.PositionLimits(2, 1) * ones(1, x_samples); % y is constant ( by plane geometry)
            z_grid = z_coord * ones(1, x_samples); % z is constant ( given parameter)

            velocity_interp = zeros(x_samples, 3);

            R_grid = [x_sample_grid; y_grid; z_grid]';

            for i = 1:x_samples % loop over all x coordinates across the channel

                Wnorm = 0; % to normalize to kernel sum when support is not full of particles

                for k = 1:obj.Npart % loop over all particles

                    r = norm(R_grid(i, :) - obj.Position(k, :));

                    if (r < 3.0 * obj.H)
                        W = obj.kernel(r);
                        Wnorm = Wnorm + W * obj.mass / obj.Density(k, 1);

                        if (obj.Type(k, 1) == 1) % fluid particle
                            velocity_interp(i, :) = velocity_interp(i, :) + obj.mass * obj.Velocity(k, :) * W / obj.Density(k, 1);
                            % else % boundary particle, velocity is stored in VelocityTransport
                            %     velocity_interp(i, :) = velocity_interp(i, :) + obj.mass * obj.VelocityTransport(k, :) * W / obj.Density(k, 1);
                        end

                    end

                end

                velocity_interp(i, :) = velocity_interp(i, :) / Wnorm;

            end

            figure(fig);
            plot(x_sample_grid, velocity_interp(:, 3), 'o-', 'DisplayName', obj.PlotName, 'MarkerSize', 4);
            xlabel('$x$'); ylabel('Velocity $z$');
            title(['SPH interpolated velocity profile at z = ' num2str(z_coord)]);
            xline(obj.PositionLimits(1, 1), 'k--', 'HandleVisibility', 'off');
            xline(obj.PositionLimits(1, 2), 'k--', 'HandleVisibility', 'off');
            xline(0.0, 'k', 'HandleVisibility', 'off');
            xline(obj.FluidBox(1), 'k', 'HandleVisibility', 'off');
            particle_initial = obj.PositionLimits(1, 1):obj.dp:obj.PositionLimits(1, 2);
            plot(particle_initial, zeros(1, length(particle_initial)), 'bo', 'HandleVisibility', 'off');

        end

        function CheckSum = CheckKernel(obj)

            CheckSum = zeros(obj.Npart, 1);

            figure; hold on;

            for k = 1:obj.Npart

                % only consider fluid particles
                if (obj.Type(k, 1) == 1)

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

        function PlotProfileOverTime(obj, fig, x_samples, z_coord, switch_var, filenums)

            for k = 1:length(filenums)
                filenum = filenums(k);
                obj = obj.Update(filenum);
                obj.PlotProfile(fig, x_samples, z_coord, switch_var);
            end

        end

        function ScatterPlotOverTime(obj, fig, filenums)

            for k = 1:length(filenums)
                filenum = filenums(k);
                obj = obj.Update(filenum);
                obj.ScatterParticlePlot(fig);
            end

        end

    end

end
