clc; clear; close all;

% Constants
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
Dt = 3;                  % Time step, minutes

% Time range
T_total = 120;  % Total simulation time
time = (0:1:T_total-1)'; % 1-min step time vector
time_steps = length(time);

% Define ambient temperatures
Tambient_values = -10:5:35; % Ambient temp range
num_cases = length(Tambient_values);

% Initialize storage
TOT_cases = NaN(time_steps, num_cases);
HST_cases = NaN(time_steps, num_cases);
final_HST = NaN(1, num_cases);
LOL_cumulative_days_all = NaN(time_steps, num_cases); % <<< Storage for cumulative LOL in days

% Load factor is constant
K = ones(time_steps, 1) * 1;

% Loop over all ambient temperatures
for case_num = 1:num_cases
    Tambient = Tambient_values(case_num); % Current ambient temp

    % Initialize
    TOT = NaN(time_steps, 1);
    HST = NaN(time_steps, 1);
    PU_life = NaN(time_steps, 1);
    LOL_cumulative_minutes = NaN(time_steps, 1);

    % Initial conditions (cold start)
    TOT(1) = Tambient;
    HST(1) = Tambient;
    delta_theta_h1 = 0;
    delta_theta_h2 = 0;

    % Iteration
    for i = 2:time_steps
        D_theta_0 = (Dt / (k11 * tao_0)) * ...
            (((1 + K(i)^2 * R) / (1 + R))^x * delta_theta_or - (TOT(i-1) - Tambient));
        TOT(i) = TOT(i-1) + D_theta_0;

        D_delta_theta_h1 = (Dt / (k22 * tao_w)) * (k21 * delta_theta_hr * K(i)^y - delta_theta_h1);
        delta_theta_h1 = delta_theta_h1 + D_delta_theta_h1;

        D_delta_theta_h2 = (Dt / ((1 / k22) * tao_0)) * ((k21 - 1) * delta_theta_hr * K(i)^y - delta_theta_h2);
        delta_theta_h2 = delta_theta_h2 + D_delta_theta_h2;

        delta_theta_h = delta_theta_h1 - delta_theta_h2;
        HST(i) = TOT(i) + delta_theta_h;
    end

    % Calculate per-unit life
    PU_life = exp(15000 / 383 - 15000 ./ (HST + 273));

    % Cumulative loss of life in minutes
    LOL_cumulative_minutes = cumsum(PU_life) * Dt;

    % Convert to days
    LOL_cumulative_days_all(:, case_num) = LOL_cumulative_minutes / (24 * 60);

    % Store final HST
    TOT_cases(:, case_num) = TOT;
    HST_cases(:, case_num) = HST;
    final_HST(case_num) = HST(end);
end

%% Plot Hot-Spot Temperatures for Different Ambient Temperatures
figure;
plot(time, HST_cases, 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)', 'FontSize', 16);
ylabel('Hot-Spot Temperature (°C)', 'FontSize', 16);
legend(arrayfun(@(T) sprintf('Tambient = %d°C', T), Tambient_values, 'UniformOutput', false), ...
    'Location', 'northwest', 'FontSize', 12);
%title('Hot-Spot Temperature vs Time for Different Ambient Temperatures', 'FontSize', 18);
set(gca, 'FontSize', 14);

%% Plot Final Hot-Spot Temperature vs Ambient Temperature
figure;
plot(Tambient_values, final_HST, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
grid on;
xlabel('Ambient Temperature (°C)', 'FontSize', 16);
ylabel('Steady State Hot-Spot Temperature (°C)', 'FontSize', 16);
%title('Final Hot-Spot Temperature vs Ambient Temperature', 'FontSize', 18);
set(gca, 'FontSize', 14);

%% Plot Cumulative Loss of Life in Days
figure;
plot(time, LOL_cumulative_days_all, 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)', 'FontSize', 16);
ylabel('Cumulative Loss of Life (days)', 'FontSize', 16);
legend(arrayfun(@(T) sprintf('Tambient = %d°C', T), Tambient_values, 'UniformOutput', false), ...
    'Location', 'northwest', 'FontSize', 12);

ylim([0 1]);             % <<< Set y-axis limit from 0 to 1
yticks(0:0.2:1);         % <<< Set y-ticks every 0.2
set(gca, 'FontSize', 14);



%title('Cumulative Loss of Life vs Time (Days)', 'FontSize', 18);
