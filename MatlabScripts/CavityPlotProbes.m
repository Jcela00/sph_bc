clear all; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);
set(groot, 'defaultFigureUnits', 'centimeters');
set(groot, 'defaultFigurePosition', [0, 0, 16.0, 16.0]); %double column

% probe 0 -> ux(y)
% probe 1 -> uy(x)
Re = [100 1000 10000];
res = [50 100 200];
nfiles = [300 500 5000];
linestyles = {':', '--', '-'};
locations = {'South', 'South', 'South'};

for n = 1:length(Re)
    current_Re = Re(n);

    t = tiledlayout(1, 1);

    % ghia data
    vxref = csvread(['../CSV_Data/Cavity/Ghia/vx_' num2str(current_Re) '.csv'], 1, 0);

    vyref = csvread(['../CSV_Data/Cavity/Ghia/vy_' num2str(current_Re) '.csv'], 1, 0);

    ax1 = axes(t);
    plot(ax1, vyref(:, 1), vyref(:, 2), 'sr', 'DisplayName', 'Ghia, $u_y$, 257x257'); hold on;
    xlabel('$x$'); ylabel('$u_y(x)$');
    ax1.Box = 'off';
    set(ax1, 'FontSize', 11);

    ax2 = axes(t);
    plot(ax2, vxref(:, 2), vxref(:, 1), 'ob', 'DisplayName', 'Ghia, $u_x$, 257x257'); hold on;
    xlabel('$u_x(y)$'); ylabel('$y$');
    ax2.Box = 'off';
    ax2.XAxisLocation = 'top';
    ax2.YAxisLocation = 'right';
    ax2.Color = 'none';
    set(ax2, 'FontSize', 11);

    if (current_Re == 10000)

        current_res = 200;

        [ux_new, y_new] = ReadProbeData(['../CSV_Data/Cavity/probes_0_Cavity_new_summation_Re' num2str(current_Re) '_' num2str(current_res) '_' num2str(current_res) '_1rf_1prc/file_' num2str(nfiles(n)) '.csv'], 2);
        [uy_new, x_new] = ReadProbeData(['../CSV_Data/Cavity/probes_1_Cavity_new_summation_Re' num2str(current_Re) '_' num2str(current_res) '_' num2str(current_res) '_1rf_1prc/file_' num2str(nfiles(n)) '.csv'], 1);

        plot(ax1, x_new, uy_new, 'k', 'LineStyle', linestyles{m}, 'DisplayName', ['SPH ' num2str(current_res)]);

        plot(ax2, ux_new, y_new, 'k', 'LineStyle', linestyles{m}, 'HandleVisibility', 'off');

    else

        for m = 1:length(res)

            current_res = res(m);

            [ux_new, y_new] = ReadProbeData(['../CSV_Data/Cavity/probes_0_Cavity_new_summation_Re' num2str(current_Re) '_' num2str(current_res) '_' num2str(current_res) '_1rf_1prc/file_' num2str(nfiles(n)) '.csv'], 2);
            [uy_new, x_new] = ReadProbeData(['../CSV_Data/Cavity/probes_1_Cavity_new_summation_Re' num2str(current_Re) '_' num2str(current_res) '_' num2str(current_res) '_1rf_1prc/file_' num2str(nfiles(n)) '.csv'], 1);

            plot(ax1, x_new, uy_new, 'k', 'LineStyle', linestyles{m}, 'DisplayName', ['SPH ' num2str(current_res)]);

            plot(ax2, ux_new, y_new, 'k', 'LineStyle', linestyles{m}, 'HandleVisibility', 'off');

        end

    end

    % Collect plot handles for both axes
    h1 = flipud(findobj(ax1.Children, '-property', 'DisplayName')); % Get all objects with a 'DisplayName' property in ax1
    axis([0.0 1.10 -0.6 0.5]);
    h2 = flipud(findobj(ax2.Children, '-property', 'DisplayName')); % Get all objects with a 'DisplayName' property in ax2
    axis([-0.4 1.0 0.0 1.0]);

    h = [h1(1, :); h2(1, :); h1(2:end, :)];

    % Combine the handles and create the legend
    l = legend([h], 'Location', locations{n});
    l.Box = 'off';
    l.FontSize = 11;

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
    exportgraphics(gcf, ['LatexFigures/CavityProbesRe' num2str(current_Re) '.pdf'], 'ContentType', 'vector', 'Resolution', 300);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [velocity, position] = ReadProbeData(filename, position_component)
    % position_component: 1 -> x, 2 -> y, 3 -> z
    probeData = csvread(filename, 1, 0);

    velocity = probeData(:, 2);
    position = probeData(:, 4:6);
    position = position(:, position_component);
    [position, idx] = sortrows(position);
    velocity = velocity(idx);
end
