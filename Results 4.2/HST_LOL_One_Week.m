clc; clear; close all;

%% Load the wind farm input data with error handling
[file, path] = uigetfile('*.csv', 'Select Wind Inputs CSV File');
if isequal(file, 0)
    error('No file selected. Please select the correct CSV file.');
else
    wind_inputs_data = readmatrix(fullfile(path, file));
end

%% Define original and new time steps
Dt_original = 10; % Original step (10 minutes)
Dt_new = 1;       % New step (1 minute)

time_original = (0:Dt_original:(length(wind_inputs_data) - 1) * Dt_original)';
time_new = (0:Dt_new:time_original(end))';

%% Perform linear interpolation
wind_inputs_interpolated = interp1(time_original, wind_inputs_data, time_new, 'linear');

%% Convert time to days for plotting
time_days_new = time_new / (24 * 60);

%% Extract 1 week (7 days) of data
one_week_minutes = 7 * 24 * 60;
indices_one_week = find(time_new <= one_week_minutes);
time_days_week = time_days_new(indices_one_week);
wind_inputs_week = wind_inputs_interpolated(indices_one_week);

% Fix invalid wind inputs
wind_inputs_week(wind_inputs_week <= 0) = 0.01;

%% Constants
delta_theta_or = 45;     % Top-oil temperature rise at rated losses, K
delta_theta_hr = 35;     % Hot-spot-to-top-oil gradient at rated current, K
tao_0 = 150;             % Average oil time constant, min
tao_w = 7;               % Winding time constant, min
R = 8;                   % Ratio of load losses at rated current to no-load losses
x = 0.8;                 % Oil exponent
y = 1.3;                 % Winding exponent
k11 = 0.5;
k21 = 2;
k22 = 2;

% Constant ambient temperature
theta_a_week = repmat(22.2, length(time_days_week), 1);

%% Initialize results
HST_week = NaN(length(time_days_week), 1);
TOT_week = NaN(length(time_days_week), 1);

% Initial conditions
K_0 = wind_inputs_week(1);
theta_0 = ((1 + K_0^2 * R) / (1 + R))^x * delta_theta_or + theta_a_week(1);
delta_theta_h1 = k21 * K_0^y * delta_theta_hr;
delta_theta_h2 = (k21 - 1) * K_0^y * delta_theta_hr;

%% Simulation loop
for i = 1:length(time_days_week)
    D_theta_0 = (Dt_new / (k11 * tao_0)) * (((1 + wind_inputs_week(i)^2 * R) / (1 + R))^x * delta_theta_or - (theta_0 - theta_a_week(i)));
    theta_0 = theta_0 + D_theta_0;

    D_delta_theta_h1 = (Dt_new / (k22 * tao_w)) * (k21 * delta_theta_hr * wind_inputs_week(i)^y - delta_theta_h1);
    delta_theta_h1 = delta_theta_h1 + D_delta_theta_h1;

    D_delta_theta_h2 = (Dt_new / ((1 / k22) * tao_0)) * ((k21 - 1) * delta_theta_hr * wind_inputs_week(i)^y - delta_theta_h2);
    delta_theta_h2 = delta_theta_h2 + D_delta_theta_h2;

    delta_theta_h = delta_theta_h1 - delta_theta_h2;
    HST_week(i) = theta_0 + delta_theta_h;
    TOT_week(i) = theta_0;
end

%% Calculate Per-Unit Life and Cumulative Loss of Life
PU_life_week = exp(15000 / 383 - 15000 ./ (HST_week + 273)); % Standard aging acceleration factor
LOL_cumulative_minutes_week = cumsum(PU_life_week) * Dt_new; % Cumulative LOL in MINUTES
LOL_cumulative_days_week = LOL_cumulative_minutes_week / (24 * 60); % Cumulative LOL in DAYS

%% Downsample for plotting (optional)
downsample_factor_week = 100;
indices_reduced_week = 1:downsample_factor_week:length(time_days_week);

%% Plot Hot-Spot and Top-Oil Temperatures
figure('Position', [100, 100, 1200, 600]);
plot(time_days_week(indices_reduced_week), HST_week(indices_reduced_week), '-r', 'LineWidth', 2); hold on;
plot(time_days_week(indices_reduced_week), TOT_week(indices_reduced_week), '-b', 'LineWidth', 2);
grid on;
xlabel('Time (Days)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 18, 'FontWeight', 'bold');
legend('Hot-Spot Temperature', 'Top-Oil Temperature', 'Location', 'best', 'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'FontSize', 18, 'FontWeight', 'bold');
hold off;

%% Plot Cumulative Loss of Life in MINUTES
figure('Position', [100, 100, 1000, 500]);
plot(time_days_week, LOL_cumulative_minutes_week, 'LineWidth', 2, 'Color', [0.2 0.5 0.8]);
grid on;
xlabel('Time (Days)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Cumulative Loss of Life (Minutes)', 'FontSize', 18, 'FontWeight', 'bold');
%title('Cumulative Loss of Life in Minutes', 'FontSize', 18, 'FontWeight', 'bold');
set(gca, 'FontSize', 16, 'FontWeight', 'bold');

% Optional limits for better visualization
% ylim([0 max(LOL_cumulative_minutes_week)*1.1]);

%% Plot Cumulative Loss of Life in DAYS
figure('Position', [100, 100, 1000, 500]);
plot(time_days_week, LOL_cumulative_days_week, 'LineWidth', 2, 'Color', [0.1 0.6 0.3]);
grid on;
xlabel('Time (Days)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Cumulative Loss of Life (Days)', 'FontSize', 16, 'FontWeight', 'bold');
%title('Cumulative Loss of Life in Days', 'FontSize', 18, 'FontWeight', 'bold');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

% Optional limits for better visualization
ylim([0 1]);
yticks(0:0.2:1);