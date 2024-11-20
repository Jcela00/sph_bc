clear all; close all; clc;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultLineLineWidth', 1.0);
set(groot, 'defaultLineMarkerSize', 2);
set(groot, 'defaultFigureUnits', 'centimeters');
% set(groot, 'defaultFigurePosition', [0, 0, 8.5, 6.0]); %single column
set(groot, 'defaultFigurePosition', [100, 100, 16.0, 10.0]); %double column

% t0 = 0.02/1.2e-4;
% t = t / t0;

paper = [0.685 0.460; 0.717 0.463];

params100 = [0.01 1 1];
params1000 = [0.001 1 1];

EllipseRe100 = ReadDragLift('../CSV_Data/dragEllipse/Re100_200.csv', 1, 1, params100);
EllipseRe100_2 = ReadDragLift('../CSV_Data/dragEllipse/Re100_400.csv', 1, 1, params100);
EllipseRe1000 = ReadDragLift('../CSV_Data/dragEllipse/Re1000_200.csv', 1, 1, params1000);
EllipseRe1000_2 = ReadDragLift('../CSV_Data/dragEllipse/Re1000_400.csv', 1, 1, params1000);

EllipseRe100{3} = EllipseRe100{3} .* EllipseRe100{2};
EllipseRe100{4} = EllipseRe100{4} .* EllipseRe100{2};
EllipseRe100_2{3} = EllipseRe100_2{3} .* EllipseRe100_2{2};
EllipseRe100_2{4} = EllipseRe100_2{4} .* EllipseRe100_2{2};

EllipseRe1000{3} = EllipseRe1000{3} .* EllipseRe1000{2};
EllipseRe1000{4} = EllipseRe1000{4} .* EllipseRe1000{2};
EllipseRe1000_2{3} = EllipseRe1000_2{3} .* EllipseRe1000_2{2};
EllipseRe1000_2{4} = EllipseRe1000_2{4} .* EllipseRe1000_2{2};

EllipseRe100{3} = LowPassFilter(EllipseRe100{3}, EllipseRe100{1}, 'Re100 Drag');
EllipseRe100{4} = LowPassFilter(EllipseRe100{4}, EllipseRe100{1}, 'Re100 Lift');
EllipseRe100_2{3} = LowPassFilter(EllipseRe100_2{3}, EllipseRe100_2{1}, 'Re100 400 Drag');
EllipseRe100_2{4} = LowPassFilter(EllipseRe100_2{4}, EllipseRe100_2{1}, 'Re100 400 Lift');

EllipseRe1000{3} = LowPassFilter(EllipseRe1000{3}, EllipseRe1000{1}, 'Re1000 Drag');
EllipseRe1000{4} = LowPassFilter(EllipseRe1000{4}, EllipseRe1000{1}, 'Re1000 Lift');
EllipseRe1000_2{3} = LowPassFilter(EllipseRe1000_2{3}, EllipseRe1000_2{1}, 'Re1000 400 Drag');
EllipseRe1000_2{4} = LowPassFilter(EllipseRe1000_2{4}, EllipseRe1000_2{1}, 'Re1000 400 Lift');

own_fc = zeros(2, 2);

% Drag plot Re 100
window = [10 38];
figure; hold on;
plot(EllipseRe100{1}, EllipseRe100{3}, 'b-', 'DisplayName', '$Re=100$ $N=50$');
plot(EllipseRe100_2{1}, EllipseRe100_2{3}, 'r-', 'DisplayName', '$Re=100$ $N=100$');
[avg] = AverageInTimeWindow(EllipseRe100{3}, EllipseRe100{1}, window);
[avg2] = AverageInTimeWindow(EllipseRe100_2{3}, EllipseRe100_2{1}, window);
fprintf('Re 100 avg Cd: %f\n', avg);
fprintf('Re 100 avg Cd: %f\n', avg2);
yline(avg, '--b', 'DisplayName', ['Avg: ' num2str(avg)]);
yline(avg2, '--r', 'DisplayName', ['Avg: ' num2str(avg2)]);
xlabel('Time [s]');
ylabel('$C_D$');
title('Drag Coefficient Re=100');
legend('Location', 'best', 'box', 'off');
set(gca, 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/EllipseDrag100.pdf'], 'ContentType', 'vector', 'Resolution', 300);

own_fc(1, 1) = avg;
own_fc(2, 1) = avg2;

% Drag plot Re 1000
figure; hold on;
plot(EllipseRe1000{1}, EllipseRe1000{3}, 'b-', 'DisplayName', '$Re=1000$ $N=50$');
plot(EllipseRe1000_2{1}, EllipseRe1000_2{3}, 'r-', 'DisplayName', '$Re=1000$ $N=100$');
[avg] = AverageInTimeWindow(EllipseRe1000{3}, EllipseRe1000{1}, window);
[avg2] = AverageInTimeWindow(EllipseRe1000_2{3}, EllipseRe1000_2{1}, window);
fprintf('Re 1000 avg Cd: %f\n', avg);
fprintf('Re 1000 avg Cd: %f\n', avg2);
yline(avg, '--b', 'DisplayName', ['Avg: ' num2str(avg)]);
yline(avg2, '--r', 'DisplayName', ['Avg: ' num2str(avg2)]);
xlabel('Time [s]');
ylabel('$C_D$');
title('Drag Coefficient Re=1000');
legend('Location', 'best', 'box', 'off');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/EllipseDrag1000.pdf'], 'ContentType', 'vector', 'Resolution', 300);

% LIFT PLOT
% Lift plot Re 100
figure; hold on;
plot(EllipseRe100{1}, EllipseRe100{4}, 'b-', 'DisplayName', '$Re=100$ N=50');
plot(EllipseRe100_2{1}, EllipseRe100_2{4}, 'r-', 'DisplayName', '$Re=100$ N=100');
[avg] = AverageInTimeWindow(EllipseRe100{4}, EllipseRe100{1}, window);
[avg2] = AverageInTimeWindow(EllipseRe100_2{4}, EllipseRe100_2{1}, window);
fprintf('Re 100 avg Cl: %f\n', avg);
fprintf('Re 100 avg Cl: %f\n', avg2);
yline(avg, '--b', 'DisplayName', ['Avg: ' num2str(avg)]);
yline(avg2, '--r', 'DisplayName', ['Avg: ' num2str(avg2)]);
xlabel('Time [s]');
ylabel('$C_L$');
title('Lift Coefficient Re=100');
legend('Location', 'best', 'box', 'off');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/EllipseLift100.pdf'], 'ContentType', 'vector', 'Resolution', 300);

own_fc(1, 2) = avg;
own_fc(2, 2) = avg2;

% Lift plot Re 1000
figure; hold on;
plot(EllipseRe1000{1}, EllipseRe1000{4}, 'b-', 'DisplayName', '$Re=1000$ N=50');
plot(EllipseRe1000_2{1}, EllipseRe1000_2{4}, 'r-', 'DisplayName', '$Re=1000$ N=100');
[avg] = AverageInTimeWindow(EllipseRe1000{4}, EllipseRe1000{1}, window);
[avg2] = AverageInTimeWindow(EllipseRe1000_2{4}, EllipseRe1000_2{1}, window);
fprintf('Re 1000 avg Cl: %f\n', avg);
fprintf('Re 1000 avg Cl: %f\n', avg2);
yline(avg, '--k', 'DisplayName', ['Avg: ' num2str(avg)]);
yline(avg2, '--r', 'DisplayName', ['Avg: ' num2str(avg2)]);
xlabel('Time [s]');
ylabel('$C_L$');
title('Lift Coefficient Re=1000');
legend('Location', 'best', 'box', 'off');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/EllipseLift1000.pdf'], 'ContentType', 'vector', 'Resolution', 300);

ratio = own_fc ./ paper;

disp(ratio);

function [DragDataset] = ReadDragLift(filename, t0, clonenormalize, params)
    data = csvread(filename, 1, 0);
    t = data(:, 1) / t0;
    u = data(:, 2);
    drag = data(:, 3);
    lift = data(:, 4);

    if (clonenormalize == 1)
        nu = params(1);
        U = params(2);
        D = params(3);
        L = 1.26315186324729;
        drag = drag * nu / (0.5 * U * D * L);
        lift = lift * nu / (0.5 * U * D * L);
    end

    DragDataset = {t, u, drag, lift};
end

function [mean_value] = AverageInTimeWindow(signal, t, window)
    mean_value = mean(signal(t > window(1) & t < window(2)));
end

function [filtered_signal] = LowPassFilter(signal, t, titlestr)
    % FILTER_NOISE Filters noise from a discrete signal in frequency space.
    %
    % INPUTS:
    % signal - The noisy input signal (vector)
    % t     - time (vector)
    %
    % OUTPUTS:
    % filtered_signal - The filtered signal (vector)
    % cutoff_frequency - The programmatically determined cutoff frequency (scalar, in Hz)

    dt = t(2:end) - t(1:end - 1);
    dt = mean(dt);
    Fs = 1 / dt; % Sampling frequency

    % Remove the DC component (mean of the signal)
    signal_detrended = signal - mean(signal);

    % Compute FFT of the signal
    L = length(signal); % Length of signal
    Y = fft(signal_detrended); % FFT of the detrended signal
    f = (0:L - 1) * (Fs / L); % Frequency vector (unshifted)

    f_cutoff = f(60);
    mask = f < f_cutoff;
    mask = mask | flip(mask);

    Y_filtered = Y .* mask';

    filtered_signal = ifft(Y_filtered, 'symmetric') + mean(signal);

    % Y_shifted = fftshift(Y); % Shift the FFT
    % f_shifted = (-L / 2:L / 2 - 1) * Fs / L; % Frequency vector (shifted)

    % % Compute the magnitude spectrum
    % magnitude_spectrum = abs(Y_shifted) .^ 2;

    % % Compute cumulative energy spectrum
    % cumulative_energy = cumsum(magnitude_spectrum);
    % total_energy = sum(magnitude_spectrum);

    % % Find the cutoff frequency corresponding to 95% of the total energy
    % cutoff_index = find(cumulative_energy >= 0.75 * total_energy, 1);
    % cutoff_frequency = f_shifted(cutoff_index);

    % % Create a low-pass filter mask
    % filter_mask = f <= cutoff_frequency; % Low-pass filter mask
    % filter_mask = filter_mask | flip(filter_mask); % Symmetry for negative frequencies

    % % Apply the filter in the frequency domain
    % Y_filtered = Y .* filter_mask';

    % % Transform back to time domain
    % filtered_signal = ifft(Y_filtered, 'symmetric') + mean(signal);

    % Display the results
    fprintf('Cutoff frequency: %.2f Hz\n', f_cutoff);

    % figure;
    % subplot(3, 1, 1);
    % plot(t, signal);
    % title(['Original Signal with Noise ' titlestr]);
    % subplot(3, 1, 2);
    % plot(t, filtered_signal);
    % title('Filtered Signal');
    % subplot(3, 1, 3);
    % semilogy(f, abs(Y) .^ 2);
    % xline(f_cutoff, '--r');
end
