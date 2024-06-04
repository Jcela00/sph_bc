classdef ParticleData

    properties
        BaseName
        PlotName
        dim
        Type
        Density
        DensityPrev
        Pressure
        DeltaDensity
        Force
        Velocity
        VelMax
        VelMin
        VelocityPrev
        Domain
        Position
        Npart
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
            obj.PlotName = [obj.BaseName ' ' num2str(filenum)];
            obj.dim = dim;
            obj.Type = data(:, 1);
            obj.Density = data(:, 2);
            obj.DensityPrev = data(:, 3);
            obj.Pressure = data(:, 4);
            obj.DeltaDensity = data(:, 5);
            obj.Force = data(:, 6:8);
            obj.Velocity = data(:, 9:11);
            obj.VelMax = max(sqrt(sum(obj.Velocity .^ 2, 2)));
            obj.VelMin = min(sqrt(sum(obj.Velocity .^ 2, 2)));
            obj.VelocityPrev = data(:, 12:14);
            obj.Domain = data(:, 15);
            obj.Position = data(:, 16:18);
            obj.Npart = length(obj.Type);
            obj.dp = dp;
            obj.dim_self = dim;
            obj.Hconst_self = Hconst;

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

        function obj = Update(obj, filenum)
            % create a new instance of the particle data object at a new time step
            obj = ParticleData(obj.filename_base, filenum, obj.BaseName, obj.dp, obj.Nfluid, obj.Nboundary, obj.rho0, obj.dim, obj.Hconst_self);
        end

        function PlotParticles(obj, fig)
            hold on;
            custom_colormap = turbo;

            for k = 1:obj.Npart
                vel_magnitude = norm(obj.Velocity(k, :));
                color_index = round(1 + ((255 - 1) / (obj.VelMax - obj.VelMin)) * (vel_magnitude - obj.VelMin));
                plot(obj.Position(k, 3), obj.Position(k, 1), 'MarkerEdgeColor', custom_colormap(color_index, :), 'MarkerFaceColor', custom_colormap(color_index, :), 'Marker', 'o', 'MarkerSize', 4);
            end

            fig.Colormap = custom_colormap;

            colorbar;
            clim([obj.VelMin, obj.VelMax]);
            axis equal;
            xlabel('$z$'); ylabel('$x$');
            title(obj.PlotName);

        end

        function ScatterParticlePlot(obj, fig)
            figure(fig);

            for k = 1:obj.Npart

                if obj.Type(k) == 1

                    if (k == 1)
                        plot(obj.Position(:, 1), obj.Velocity(:, 3), 'o', 'MarkerSize', 2, 'DisplayName', obj.PlotName);
                    else
                        plot(obj.Position(:, 1), obj.Velocity(:, 3), 'o', 'MarkerSize', 2, 'HandleVisibility', 'off');
                    end

                else

                    if (k == 1)
                        plot(obj.Position(:, 1), obj.VelocityPrev(:, 3), 'o', 'MarkerSize', 2, 'DisplayName', obj.PlotName);

                    else
                        plot(obj.Position(:, 1), obj.VelocityPrev(:, 3), 'o', 'MarkerSize', 2, 'HandleVisibility', 'off');
                    end

                end

            end

            xlabel('$x$'); ylabel('Velocity $z$');
            xline(obj.PositionLimits(1, 1), 'k', 'HandleVisibility', 'off');
            xline(obj.PositionLimits(1, 2), 'k--', 'HandleVisibility', 'off');
            xline(0.0, 'k', 'HandleVisibility', 'off');
            xline(obj.FluidBox(1), 'k', 'HandleVisibility', 'off');
            particle_initial = obj.PositionLimits(1, 1):obj.dp:obj.PositionLimits(1, 2);
            plot(particle_initial, zeros(1, length(particle_initial)), 'bo', 'HandleVisibility', 'off');
            axis equal;

        end

        function W = kernel(obj, r)
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

        function PlotProfile(obj, fig, x_samples, z_coord, switch_var)

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

                    if (r < 2 * obj.H)
                        W = obj.kernel(r);
                        Wnorm = Wnorm + W * obj.mass / obj.Density(k, 1);

                        if (obj.Type(k, 1) == 1) % fluid particle
                            velocity_interp(i, :) = velocity_interp(i, :) + obj.mass * obj.Velocity(k, :) * W / obj.Density(k, 1);
                        else % boundary particle, velocity is stored in VelocityPrev
                            velocity_interp(i, :) = velocity_interp(i, :) + obj.mass * obj.VelocityPrev(k, :) * W / obj.Density(k, 1);
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
            axis equal;

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
                filenum = filenums(k)
                obj = obj.Update(filenum);
                obj.ScatterParticlePlot(fig);
            end

        end

    end

end
