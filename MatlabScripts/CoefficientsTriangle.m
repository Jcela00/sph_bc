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

TriangleRe100_25 = ReadDragLift('../CSV_Data/DragTriangle/Re100_25.csv', 1, 1, params100);
TriangleRe100_50 = ReadDragLift('../CSV_Data/DragTriangle/Re100_50.csv', 1, 1, params100);

TriangleRe1000_25 = ReadDragLift('../CSV_Data/DragTriangle/Re1000_25.csv', 1, 1, params1000);
TriangleRe1000_50 = ReadDragLift('../CSV_Data/DragTriangle/Re1000_50.csv', 1, 1, params1000);

%cols 1 2 3 4 | drag average | drag amplitude | drag frequency1 | drag frequency2 | lift average | lift amplitude | lift frequency1 | lift frequency2
Re100_table = zeros(2, 8);
Re1000_table = zeros(2, 8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RE 100 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% window = [55.9 95.47];
% fig1 = figure; hold on;
% Re100_table = PlotCoefficient(TriangleRe100_25, window, fig1, Re100_table, 1, 'b', 0, 'Re100_25');
% % exportgraphics(gcf, ['LatexFigures/TriangleDrag100_25.pdf'], 'ContentType', 'vector', 'Resolution', 300);

window = [58.5 90.88];
fig2 = figure; hold on;
Re100_table = PlotCoefficient(TriangleRe100_50, window, fig2, Re100_table, 2, 'r', 0, 'Re100_50');
exportgraphics(gcf, ['LatexFigures/TriangleDrag100_50.pdf'], 'ContentType', 'vector', 'Resolution', 300);

% window = [56.2 93.4];
% fig3 = figure; hold on;
% Re100_table = PlotCoefficient(TriangleRe100_25, window, fig3, Re100_table, 1, 'b', 1, 'Re100_25');
% % exportgraphics(gcf, ['LatexFigures/TriangleLift100_25.pdf'], 'ContentType', 'vector', 'Resolution', 300);

window = [58.95 91.3];
fig4 = figure; hold on;
Re100_table = PlotCoefficient(TriangleRe100_50, window, fig4, Re100_table, 2, 'r', 1, 'Re100_50');
exportgraphics(gcf, ['LatexFigures/TriangleLift100_50.pdf'], 'ContentType', 'vector', 'Resolution', 300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RE 1000 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% window = [55.9 95.47];
% fig5 = figure; hold on;
% Re1000_table = PlotCoefficient(TriangleRe1000_25, window, fig5, Re1000_table, 1, 'b', 0, 'Re1000_25');
% % exportgraphics(gcf, ['LatexFigures/TriangleDrag1000_25.pdf'], 'ContentType', 'vector', 'Resolution', 300);

window = [58.5 90.88];
fig6 = figure; hold on;
Re1000_table = PlotCoefficient(TriangleRe1000_50, window, fig6, Re1000_table, 2, 'r', 0, 'Re1000_50');
exportgraphics(gcf, ['LatexFigures/TriangleDrag1000_50.pdf'], 'ContentType', 'vector', 'Resolution', 300);

% window = [50.2 91.8];
% fig7 = figure; hold on;
% Re1000_table = PlotCoefficient(TriangleRe1000_25, window, fig7, Re1000_table, 1, 'b', 1, 'Re1000_25');
% % exportgraphics(gcf, ['LatexFigures/TriangleLift1000_25.pdf'], 'ContentType', 'vector', 'Resolution', 300);

window = [49.8 91.8];
fig8 = figure; hold on;
Re1000_table = PlotCoefficient(TriangleRe1000_50, window, fig8, Re1000_table, 2, 'r', 1, 'Re1000_50');
exportgraphics(gcf, ['LatexFigures/TriangleLift1000_50.pdf'], 'ContentType', 'vector', 'Resolution', 300);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TableRe100 = table(Re100_table(:, 1), Re100_table(:, 2), Re100_table(:, 3), Re100_table(:, 4), Re100_table(:, 5), Re100_table(:, 6), Re100_table(:, 7), Re100_table(:, 8));
TableRe100.Properties.VariableNames = ["DragAvg", "DragAmp", "DragFreq1", "DragFreq2", "LiftAvg", "LiftAmp", "LiftFreq1", "LiftFreq2"];
TableRe100.Properties.RowNames = {'N=25', 'N=50'};

TableRe1000 = table(Re1000_table(:, 1), Re1000_table(:, 2), Re1000_table(:, 3), Re1000_table(:, 4), Re1000_table(:, 5), Re1000_table(:, 6), Re1000_table(:, 7), Re1000_table(:, 8));
TableRe1000.Properties.VariableNames = ["DragAvg", "DragAmp", "DragFreq1", "DragFreq2", "LiftAvg", "LiftAmp", "LiftFreq1", "LiftFreq2"];
TableRe1000.Properties.RowNames = {'N=25', 'N=50'};

Re100_table_paper = [1.95 0.075 0.214 * 2 0.214 * 2 0.0 0.473 0.214 0.214];
TableRe100Paper = table(Re100_table_paper(:, 1), Re100_table_paper(:, 2), Re100_table_paper(:, 3), Re100_table_paper(:, 4), Re100_table_paper(:, 5), Re100_table_paper(:, 6), Re100_table_paper(:, 7), Re100_table_paper(:, 8));
TableRe100Paper.Properties.VariableNames = ["DragAvg", "DragAmp", "DragFreq1", "DragFreq2", "LiftAvg", "LiftAmp", "LiftFreq1", "LiftFreq2"];
TableRe100Paper.Properties.RowNames = {'Exact'};

Re1000_table_paper = [2.33 0.33 0.213 * 2 0.213 * 2 0.0 0.65 0.213 0.213];

TableRe1000Paper = table(Re1000_table_paper(:, 1), Re1000_table_paper(:, 2), Re1000_table_paper(:, 3), Re1000_table_paper(:, 4), Re1000_table_paper(:, 5), Re1000_table_paper(:, 6), Re1000_table_paper(:, 7), Re1000_table_paper(:, 8));
TableRe1000Paper.Properties.VariableNames = ["DragAvg", "DragAmp", "DragFreq1", "DragFreq2", "LiftAvg", "LiftAmp", "LiftFreq1", "LiftFreq2"];
TableRe1000Paper.Properties.RowNames = {'Exact'};

disp('-----------------------------------------------');
disp('Computed values for Re100');
disp(TableRe100);
disp('Computed values for Re1000');
disp(TableRe1000);
disp('-----------------------------------------------');

disp('Paper values for Re100');
disp(TableRe100Paper);
disp('Paper values for Re1000');
disp(TableRe1000Paper);

ratio100 = Re100_table ./ [Re100_table_paper; Re100_table_paper];
ratio1000 = Re1000_table ./ [Re1000_table_paper; Re1000_table_paper];

disp('Ratio values for Re100');
disp(ratio100);
disp('Ratio values for Re1000');
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
        L = 1.0;
        drag = drag / (0.5 * U * D * L);
        lift = lift / (0.5 * U * D * L);
    end

    % u = u(t > window(1) & t < window(2));
    % drag = drag(t > window(1) & t < window(2));
    % lift = lift(t > window(1) & t < window(2));
    % t = t(t > window(1) & t < window(2));

    [drag] = LowPassFilter(drag, t);
    [lift] = LowPassFilter(lift, t);

    DragDataset = {t, u, drag, lift};

end

function [mean_value, amplitude, frequencies] = AverageInTimeWindow(signal, t, window)

    % compute mean value
    mean_value = mean(signal(t > window(1) & t < window(2)));

    % compute amplitude
    t_cut = t(t > window(1) & t < window(2));

    dt = t_cut(2:end) - t_cut(1:end - 1);
    dt = mean(dt);

    signal_cut = signal(t > window(1) & t < window(2));
    signal_cut = signal_cut - mean(signal_cut);
    % amplitude = (max(signal_cut) - min(signal_cut)) / 2;

    % compute amplitude as average of multiple maximum - minimum / 2
    % find peaks of signal
    peaks_max_idx = find(signal_cut(2:end - 1) > signal_cut(1:end - 2) & signal_cut(2:end - 1) > signal_cut(3:end)) + 1;
    peaks_min_idx = find(signal_cut(2:end - 1) < signal_cut(1:end - 2) & signal_cut(2:end - 1) < signal_cut(3:end)) + 1;

    maxvalues = signal_cut(peaks_max_idx);
    minvalues = signal_cut(peaks_min_idx);

    % some values are just local extrema (this is more true for very noise signals), we take the first 10 as real extrema
    lastidx = min(min(length(maxvalues), length(minvalues)), 6);
    maxvalues = sort(maxvalues); maxvalues = flip(maxvalues); maxvalues = maxvalues(1:lastidx);
    minvalues = sort(minvalues); minvalues = minvalues(1:lastidx);

    amplitude = (mean(maxvalues) - mean(minvalues)) / 2;

    % compute amplitude through rms

    amplitude_rms = sqrt(2) * sqrt(mean(signal_cut .^ 2));
    amplitude = amplitude_rms;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZE EXTREMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % max_idx = find(ismember(signal_cut, maxvalues));
    % min_idx = find(ismember(signal_cut, minvalues));
    % figure;
    % plot(t_cut, signal_cut); hold on;
    % plot(t_cut(max_idx), signal_cut(max_idx), 'or', 'MarkerSize', 10);
    % plot(t_cut(min_idx), signal_cut(min_idx), 'ob', 'MarkerSize', 10);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZE EXTREMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute frequency through autocorrelation
    [R, lags] = xcorr(signal_cut, 'coeff');
    lags = lags * dt; % Convert lags to time

    %%%%%% UNCOMMENT TO PLOT AUTOCORRELATION %%%%%%%%%%%%%%%%%%%%%
    % fig1 = figure;
    % plot(lags, R);
    % xlabel('Lag [s]');
    % ylabel('Autocorrelation');
    % pause
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % find peaks
    R_positive = R(lags > 0);
    lags_positive = lags(lags > 0);

    peaks_idx = find(R_positive(2:end - 1) > R_positive(1:end - 2) & R_positive(2:end - 1) > R_positive(3:end)) + 1;

    %%%%%%%%%%%%%%%%%%% UNCOMMENT TO PLOT PEAKS %%%%%%%%%%%%%%%%%%%%%
    % figure;
    % plot(lags_positive, R_positive);
    % hold on;
    % plot(lags_positive(peaks_idx), R_positive(peaks_idx), 'ro');
    % xlabel('Lag [s]');
    % ylabel('Autocorrelation');
    % pause
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isempty(peaks_idx)
        T_dominant = lags_positive(peaks_idx(1));
        f_xcorr = 1 / T_dominant;
    else
        T_dominant = 0;
        f_xcorr = 0;
    end

    % compute frequency through zero crossings
    zero_crossings = find(diff(sign(signal_cut)) ~= 0);
    time_between_zero_crossings = t_cut(zero_crossings);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNCOMMENT TO PLOT ZERO CROSSINGS %%%%%%%%%%%%%%%%%%%%%
    % figure;
    % plot(t_cut, signal_cut);
    % hold on;
    % plot(time_between_zero_crossings, zeros(size(time_between_zero_crossings)), 'ro');
    % xlabel('Time [s]');
    % ylabel('Signal');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % compute periods
    T = diff(time_between_zero_crossings);
    T_mean = 2 * mean(T);
    f_zerocrossing = 1 / T_mean;

    % fprintf('Frequency through zero crossings: %.2f Hz, Period %.2f s\n', f_dominant, T_mean);
    % fprintf('Frequency through autocorrelation: %.2f Hz, Period %.2f s\n', f_dominant_xcorr, T_dominant);
    frequencies = [f_zerocrossing, f_xcorr];
    % frequencies = [amplitude_rms, f_xcorr];

end

function [filtered_signal] = LowPassFilter(signal, t)
    % window signal

    dt = t(2:end) - t(1:end - 1);
    dt = mean(dt);
    Fs = 1 / dt;

    L = length(signal);
    Y = fft(signal);
    f = [0:L - 1] * (Fs / L);

    f_cutoff = f(200);
    mask = f < f_cutoff;
    mask = mask | flip(mask);
    Yfiltered = Y .* mask';
    filtered_signal = ifft(Yfiltered, 'symmetric');

    % select 100 more relevant modes
    % [~, idx] = maxk(abs(Y), 50);
    % mask = ismember(1:L, idx);
    % mask = mask;
    % Yfiltered = Y .* mask';
    % filtered_signal = ifft(Yfiltered, 'symmetric');

end

function [return_table] = PlotCoefficient(dataset, window, fighandle, coefficient_table, tablerow, clr, DragOrLift, name)

    % extract reynolds number and resolution from name
    % name is formated as ReXXX_YY where XXX is the Reynolds number and YY is the resolution in unitos of D/Delta x
    splitName = split(name, '_');
    Re = str2double(splitName{1}(3:end));
    N = str2double(splitName{2});

    txt = {['$Re = ' num2str(Re) '$'], ['$N = ' num2str(N) '$']};

    if DragOrLift == 0
        index = 3;
        ylbl = '$C_d$';
        ttl = 'Drag Coefficient Re = ';
        cols = [1 4];
    elseif DragOrLift == 1
        index = 4;
        ylbl = '$C_l$';
        ttl = 'Lift Coefficient Re = ';
        cols = [5 8];
    end

    figure(fighandle)

    plot(dataset{1}, dataset{index}, '-', 'Color', clr, 'DisplayName', ['$N=' num2str(N) '$']);
    [avg, amp, freqs] = AverageInTimeWindow(dataset{index}, dataset{1}, window);

    figure(fighandle);
    % yline(avg, '--', 'Color', clr, 'DisplayName', ['Avg: ' num2str(avg) ' Amplitude: ' num2str(amp)]);
    % yline(avg + amp, '--r', 'HandleVisibility', 'off');
    % yline(avg - amp, '--r', 'HandleVisibility', 'off');
    % xline(window(1), '--k', 'HandleVisibility', 'off');
    % xline(window(2), '--k', 'HandleVisibility', 'off');
    xlabel('Time [s]');
    ylabel(ylbl);
    % title([ttl num2str(Re)])
    % legend('Location', 'best', 'box', 'off');
    set(gca, 'FontSize', 11);

    coefficient_table(tablerow, cols(1):cols(2)) = [avg amp freqs];

    return_table = coefficient_table;

end
