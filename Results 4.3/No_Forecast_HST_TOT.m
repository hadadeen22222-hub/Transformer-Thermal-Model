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

% === STEP 3: Limit to First 7 Days ONLY ===
max_minutes = 7 * 1440;  % 10080 minutes = 168 hours
if length(wind_inputs_interpolated) >= max_minutes
    K_all_7days = wind_inputs_interpolated(1:max_minutes);
else
    error('Input data is shorter than 7 days!');
end
Tamb = 22.2;
N = length(K_all_7days);

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

% === STEP 5: Run Thermal Simulation for 7 Days ONLY ===
theta_0 = ((1 + K_all_7days(1)^2 * R)/(1 + R))^x * delta_theta_or + Tamb;
delta_theta_h1 = k21 * K_all_7days(1)^y * delta_theta_hr;
delta_theta_h2 = (k21 - 1) * K_all_7days(1)^y * delta_theta_hr;

HST_7days = NaN(N, 1);
TOT_7days = NaN(N, 1);

for i = 1:N
    D_theta_0 = (Dt_new / (k11 * tao_0)) * (((1 + K_all_7days(i)^2 * R)/(1 + R))^x * delta_theta_or - (theta_0 - Tamb));
    theta_0 = theta_0 + D_theta_0;

    D_delta_theta_h1 = (Dt_new / (k22 * tao_w)) * (k21 * delta_theta_hr * K_all_7days(i)^y - delta_theta_h1);
    delta_theta_h1 = delta_theta_h1 + D_delta_theta_h1;

    D_delta_theta_h2 = (Dt_new / ((1 / k22) * tao_0)) * ((k21 - 1) * delta_theta_hr * K_all_7days(i)^y - delta_theta_h2);
    delta_theta_h2 = delta_theta_h2 + D_delta_theta_h2;

    delta_theta_h = delta_theta_h1 - delta_theta_h2;
    HST_7days(i) = theta_0 + delta_theta_h;
    TOT_7days(i) = theta_0;
end

% === STEP 6: Calculate Loss of Life ===
PU_life_7days = exp(15000 / 383 - 15000 ./ (HST_7days + 273));
LOL_minutes_7days = cumsum(PU_life_7days * Dt_new);
LOL_days_7days = LOL_minutes_7days / (24 * 60);

% === STEP 7: Time Vector for Plotting ===
time_min_7days = (0:N-1);
time_hours_7days = time_min_7days / 60;

% === STEP 8: Plot HST and TOT with Fixed Threshold Zones ===

% --- Set fixed thresholds ---
threshold_safe = 100;   % Start of caution zone
threshold_warn = 120;   % Start of danger zone

% --- Set y-limits to fully include all thresholds and values ---
y_min = min([HST_7days; TOT_7days], [], 'omitnan');
y_max = max([HST_7days; TOT_7days; threshold_warn + 5]);  % small margin above danger

% --- Plot ---
figure('Name', 'HST and TOT with Fixed Zones', 'Color', 'w'); hold on;

% Safe zone (green)
fill([time_hours_7days(1), time_hours_7days(end), time_hours_7days(end), time_hours_7days(1)], ...
     [y_min, y_min, threshold_safe, threshold_safe], ...
     [0.8, 1, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Caution zone (yellow)
fill([time_hours_7days(1), time_hours_7days(end), time_hours_7days(end), time_hours_7days(1)], ...
     [threshold_safe, threshold_safe, threshold_warn, threshold_warn], ...
     [1, 1, 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Danger zone (red)
fill([time_hours_7days(1), time_hours_7days(end), time_hours_7days(end), time_hours_7days(1)], ...
     [threshold_warn, threshold_warn, y_max, y_max], ...
     [1, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plot HST and TOT
plot(time_hours_7days, HST_7days, 'r-', 'LineWidth', 1.5);
plot(time_hours_7days, TOT_7days, 'b--', 'LineWidth', 1.2);

% Final plot settings
legend('Safe Zone', 'Caution Zone', 'Danger Zone', 'HST', 'TOT', 'Location', 'best');
ylim([y_min, y_max]);
xlim([0, 168]);
grid on;

% === Make tick numbers and labels bigger ===
set(gca, 'FontSize', 15); % Bigger numbers (tick labels)
xlabel('Time (Hours)', 'FontSize', 17, 'FontWeight'); % Bigger x-axis label
ylabel('Temperature (°C)', 'FontSize', 17, 'FontWeight'); % Bigger y-axis label