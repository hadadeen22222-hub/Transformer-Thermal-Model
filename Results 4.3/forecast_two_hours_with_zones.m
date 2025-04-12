clc; clear; close all;

% === STEP 1: Load Wind Input Data ===
[file, path] = uigetfile('*.csv', 'Select Wind Inputs CSV File');
if isequal(file, 0)
    error('No file selected.');
else
    wind_inputs_data = readmatrix(fullfile(path, file));
end

% === STEP 2: Preprocess and Interpolate to 1-Minute Resolution ===
Dt_original = 10;
Dt_new = 1;
time_original = (0:Dt_original:(length(wind_inputs_data)-1)*Dt_original)';
time_new = (0:Dt_new:time_original(end))';
wind_inputs_interpolated = interp1(time_original, wind_inputs_data, time_new, 'linear');
wind_inputs_interpolated = wind_inputs_interpolated * 2;
wind_inputs_interpolated(wind_inputs_interpolated <= 0) = 0.01;

% === STEP 3: Limit to First 7 Days (real) + 2-hour forecast ===
total_minutes = 170 * 60;        % Up to 170 hours
cutoff_minutes = 168 * 60;       % Forecast starts at 168h

if length(wind_inputs_interpolated) >= cutoff_minutes
    wind_inputs_real = wind_inputs_interpolated(1:cutoff_minutes);
else
    error('Insufficient input data.');
end

% === STEP 4: Transformer Thermal Model Constants ===
delta_theta_or = 45;
delta_theta_hr = 35;
tao_0 = 150;
tao_w = 7;
R = 8;
x = 0.8;
y = 1.3;
k11 = 0.5;
k21 = 2;
k22 = 2;
Tamb = 22.2;

% === STEP 5: Sweep K values for forecast window (168h–170h) ===
K_test_values = 0.1:0.2:1.7;
num_K = length(K_test_values);
forecast_minutes = 2 * 60; % 2 hours
N = total_minutes;

HST_cases = NaN(N, num_K);
TOT_cases = NaN(N, num_K);
time_vector = (0:N-1)' / 60;  % in hours

for k_idx = 1:num_K
    % Construct full K vector: real + constant forecast
    K_const = K_test_values(k_idx);
    K_full = [wind_inputs_real; K_const * ones(forecast_minutes, 1)];

    % Initialize
    theta_0 = ((1 + K_full(1)^2 * R)/(1 + R))^x * delta_theta_or + Tamb;
    delta_theta_h1 = k21 * K_full(1)^y * delta_theta_hr;
    delta_theta_h2 = (k21 - 1) * K_full(1)^y * delta_theta_hr;

    HST = NaN(N,1);
    TOT = NaN(N,1);

    for i = 1:N
        D_theta_0 = (Dt_new / (k11 * tao_0)) * (((1 + K_full(i)^2 * R)/(1 + R))^x * delta_theta_or - (theta_0 - Tamb));
        theta_0 = theta_0 + D_theta_0;

        D_delta_theta_h1 = (Dt_new / (k22 * tao_w)) * (k21 * delta_theta_hr * K_full(i)^y - delta_theta_h1);
        delta_theta_h1 = delta_theta_h1 + D_delta_theta_h1;

        D_delta_theta_h2 = (Dt_new / ((1 / k22) * tao_0)) * ((k21 - 1) * delta_theta_hr * K_full(i)^y - delta_theta_h2);
        delta_theta_h2 = delta_theta_h2 + D_delta_theta_h2;

        delta_theta_h = delta_theta_h1 - delta_theta_h2;
        HST(i) = theta_0 + delta_theta_h;
        TOT(i) = theta_0;
    end

    HST_cases(:, k_idx) = HST;
    TOT_cases(:, k_idx) = TOT;
end

% === STEP 6: Plot Enhanced HST with Fixed Threshold Zones ===

% --- Set fixed thresholds ---
threshold_safe = 100;  % Start of caution zone
threshold_warn = 120;  % Start of danger zone

% --- Determine y-axis limits ---
all_HST_vals = HST_cases(:);
y_min = min(all_HST_vals, [], 'omitnan');
y_max = max([all_HST_vals; threshold_warn + 5]);  % Add margin above danger

% --- Create plot ---
figure('Position', [100, 100, 900, 500]); hold on;

% --- Fill Threshold Zones between 166h and 172h ---
% Safe zone (green)
fill([166 172 172 166], ...
     [y_min y_min threshold_safe threshold_safe], ...
     [0.8 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');

% Caution zone (yellow)
fill([166 172 172 166], ...
     [threshold_safe threshold_safe threshold_warn threshold_warn], ...
     [1 1 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');

% Danger zone (red)
fill([166 172 172 166], ...
     [threshold_warn threshold_warn y_max y_max], ...
     [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');

% --- Plot HST curves ---
plot(time_vector, HST_cases, 'LineWidth', 2);

% --- Highlight forecast window (168h–170h) ---
y_limits = ylim;
fill([168 170 170 168], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
     [0.9 0.9 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');

% --- Add vertical lines at forecast boundaries ---
xline(168, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(170, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% --- Labels, legends, and final settings ---
xlabel('Time (Hours)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Hot-Spot Temperature (°C)', 'FontSize', 18, 'FontWeight', 'bold');
legend(arrayfun(@(k) sprintf('K = %.1f', k), K_test_values, 'UniformOutput', false), ...
       'Location', 'northwest', 'FontSize', 16); % bigger legend font

set(gca, 'FontSize', 16); % Bigger tick numbers
xlim([166 172]);
ylim([y_min y_max]);
grid on;