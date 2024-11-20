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

params100 = [0.01 1 1];
params1000 = [0.001 1 1];

TriangleRe100 = ReadDragLift('../CSV_Data/dragTriangle/Re100_200.csv', 1, 1, params100);
TriangleRe1000 = ReadDragLift('../CSV_Data/dragTriangle/Re1000_200.csv', 1, 1, params1000);

TriangleRe100{3} = TriangleRe100{3} .* TriangleRe100{2};
TriangleRe100{4} = TriangleRe100{4} .* TriangleRe100{2};

TriangleRe1000{3} = TriangleRe1000{3} .* TriangleRe1000{2};
TriangleRe1000{4} = TriangleRe1000{4} .* TriangleRe1000{2};

TriangleRe100{3} = LowPassFilter(TriangleRe100{3}, TriangleRe100{1}, 'Re100 Drag');
TriangleRe100{4} = LowPassFilter(TriangleRe100{4}, TriangleRe100{1}, 'Re100 Lift');

TriangleRe1000{3} = LowPassFilter(TriangleRe1000{3}, TriangleRe1000{1}, 'Re1000 Drag');
TriangleRe1000{4} = LowPassFilter(TriangleRe1000{4}, TriangleRe1000{1}, 'Re1000 Lift');

own_vals_100 = zeros(1, 2);
own_vals_1000 = zeros(1, 2);

% Drag plot Re 100
window = [15 35];
figure; hold on;
plot(TriangleRe100{1}, TriangleRe100{3}, 'b-', 'DisplayName', '$Re=100$ $N=50$');
[avg, amplitude] = AverageInTimeWindow(TriangleRe100{3}, TriangleRe100{1}, window);
fprintf('Re 100 avg Cd: %f Amplitude %f\n', avg, amplitude);
yline(avg, '--b', 'DisplayName', ['Avg: ' num2str(avg)]);
yline(avg + amplitude, '--r', 'DisplayName', ['Amplitude: ' num2str(amplitude)]);
yline(avg - amplitude, '--r');
xlabel('Time [s]');
ylabel('$C_D$');
title('Drag Coefficient Re=100');
legend('Location', 'best', 'box', 'off');
set(gca, 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/TriangleDrag100.pdf'], 'ContentType', 'vector', 'Resolution', 300);
own_vals_100(1) = avg;

% Drag plot Re 1000
figure; hold on;
plot(TriangleRe1000{1}, TriangleRe1000{3}, 'b-', 'DisplayName', '$Re=1000$ $N=50$');
[avg, amplitude] = AverageInTimeWindow(TriangleRe1000{3}, TriangleRe1000{1}, window);
fprintf('Re 1000 avg Cd: %f Amplitude %f\n', avg, amplitude);
yline(avg, '--b', 'DisplayName', ['Avg: ' num2str(avg)]);
yline(avg + amplitude, '--r', 'DisplayName', ['Amplitude: ' num2str(amplitude)]);
yline(avg - amplitude, '--r');
xlabel('Time [s]');
ylabel('$C_D$');
title('Drag Coefficient Re=1000');
legend('Location', 'best', 'box', 'off');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/TriangleDrag1000.pdf'], 'ContentType', 'vector', 'Resolution', 300);
own_vals_1000(1) = avg;

% LIFT PLOT
% Lift plot Re 100
figure; hold on;
plot(TriangleRe100{1}, TriangleRe100{4}, 'b-', 'DisplayName', '$Re=100$ N=50');
[avg, amplitude] = AverageInTimeWindow(TriangleRe100{4}, TriangleRe100{1}, window);
fprintf('Re 100 avg Cl: %f Amplitude %f\n', avg, amplitude);
yline(avg, '--b', 'DisplayName', ['Avg: ' num2str(avg)]);
yline(avg + amplitude, '--r', 'DisplayName', ['Amplitude: ' num2str(amplitude)]);
yline(avg - amplitude, '--r');
xlabel('Time [s]');
ylabel('$C_L$');
title('Lift Coefficient Re=100');
legend('Location', 'best', 'box', 'off');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/TriangleLift100.pdf'], 'ContentType', 'vector', 'Resolution', 300);
own_vals_100(2) = avg;

% Lift plot Re 1000
figure; hold on;
plot(TriangleRe1000{1}, TriangleRe1000{4}, 'b-', 'DisplayName', '$Re=1000$ N=50');
[avg, amplitude] = AverageInTimeWindow(TriangleRe1000{4}, TriangleRe1000{1}, window);
fprintf('Re 1000 avg Cl: %f Amplitude %f\n', avg, amplitude);
yline(avg, '--k', 'DisplayName', ['Avg: ' num2str(avg)]);
yline(avg + amplitude, '--r', 'DisplayName', ['Amplitude: ' num2str(amplitude)]);
yline(avg - amplitude, '--r');
xlabel('Time [s]');
ylabel('$C_L$');
title('Lift Coefficient Re=1000');
legend('Location', 'best', 'box', 'off');
set(gca, 'FontSize', 11);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
exportgraphics(gcf, ['LatexFigures/TriangleLift1000.pdf'], 'ContentType', 'vector', 'Resolution', 300);
own_vals_1000(2) = avg;

paper_vals_100 = [1.92];
paper_vals_1000 = [2.87];

ratio100 = own_vals_100(1) ./ paper_vals_100;
ratio1000 = own_vals_1000(1) ./ paper_vals_1000;

disp(ratio100);
disp(ratio1000);

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
        L = 1.2;
        drag = drag * nu / (0.5 * U * D * L);
        lift = lift * nu / (0.5 * U * D * L);
    end

    DragDataset = {t, u, drag, lift};
end

function [mean_value, amplitude] = AverageInTimeWindow(signal, t, window)
    mean_value = mean(signal(t > window(1) & t < window(2)));

    signal_cut = signal(t > window(1) & t < window(2));
    signal_cut = signal_cut - mean(signal_cut);
    amplitude = (max(signal_cut) - min(signal_cut)) / 2;

end

% function [filtered_signal] = LowPassFilter(signal, t, titlestr)

%     dt = t(2:end) - t(1:end - 1);
%     dt = mean(dt);

%     fs = 1 / dt;

%     y_detrended = signal - mean(signal);

%     % Perform FFT
%     Y = fft(y_detrended);
%     N = length(signal);
%     f = (0:N - 1) * (fs / N); % Frequency vector

%     f_shifted = [-N / 2:N / 2 - 1] * fs / N;
%     Y_shifted = fftshift(Y);

%     energy = abs(Y_shifted) .^ 2;
%     cumulative_energy = cumsum(energy);
%     total_energy = sum(energy);
%     plot(f_shifted, cumulative_energy / total_energy);
%     pause
%     cutoff_energy = 0.8 * total_energy;
%     cutoff_index = find(cumulative_energy > cutoff_energy, 1);
%     f_cutoff = f_shifted(cutoff_index);

%     fprintf('Cutoff frequency: %f\n', f_cutoff);

%     % filter_mask = f < f_cutoff;
%     % filter_mask = filter_mask | flip(filter_mask);

%     filter_mask = abs(f_shifted) < f_cutoff;

%     % Apply the mask and perform inverse FFT
%     Y_filtered = Y_shifted .* filter_mask';
%     filtered_signal = ifft(ifftshift(Y_filtered), 'symmetric') + mean(signal);

%     % Plot
%     figure;
%     subplot(3, 1, 1);
%     plot(t, signal);
%     title(['Original Signal with Noise ' titlestr]);
%     subplot(3, 1, 2);
%     plot(t, filtered_signal);
%     title('Filtered Signal');
%     subplot(3, 1, 3);
%     semilogy(f_shifted, abs(fftshift(Y)) .^ 2);
%     xline(f_cutoff, '--r');
%     xline(-f_cutoff, '--r');

% end

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
