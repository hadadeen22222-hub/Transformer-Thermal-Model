clc; clear; close all;

% Constants
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
Tambient = 20;
Dt = 3; % time step in minutes

% Time setup
time_steps = 41;
time = (0:Dt:(time_steps - 1) * Dt)';
theta_a = ones(time_steps, 1) * 30; % ambient temperature profile

% K values to evaluate
K_values = 0.1:0.2:1.7;
num_cases = length(K_values);

% Storage for all cumulative LOL
LOL_cumulative_days_all = NaN(time_steps, num_cases);
LOL_cumulative_minutes_all = NaN(time_steps, num_cases);

% Loop over K values
for case_num = 1:num_cases
    K = ones(time_steps, 1) * K_values(case_num);

    % Initialize
    TOT = NaN(time_steps, 1);
    HST = NaN(time_steps, 1);
    PU_life = NaN(time_steps, 1);
    LOL_cumulative_hours = NaN(time_steps, 1); % internal cumulative in hours

    % Initial conditions
    theta_0 = 20;
    delta_theta_h1 = 0;
    delta_theta_h2 = 0;

    for i = 1:time_steps
        % Top-oil temperature
        D_theta_0 = (Dt / (k11 * tao_0)) * ...
            (((1 + K(i)^2 * R) / (1 + R))^x * delta_theta_or - (theta_0 - theta_a(i)));
        theta_0 = theta_0 + D_theta_0;

        % Hot-spot gradient components
        D_delta_theta_h1 = (Dt / (k22 * tao_w)) * ...
            (k21 * delta_theta_hr * K(i)^y - delta_theta_h1);
        delta_theta_h1 = delta_theta_h1 + D_delta_theta_h1;

        D_delta_theta_h2 = (Dt / ((1 / k22) * tao_0)) * ...
            ((k21 - 1) * delta_theta_hr * K(i)^y - delta_theta_h2);
        delta_theta_h2 = delta_theta_h2 + D_delta_theta_h2;

        delta_theta_h = delta_theta_h1 - delta_theta_h2;

        TOT(i) = theta_0;
        HST(i) = theta_0 + delta_theta_h;

        PU_life(i) = exp(15000 / 383 - 15000 / (HST(i) + 273)); % per-unit life consumption

        % Cumulative loss of life in hours
        LOL_cumulative_hours(i) = sum(PU_life(1:i) * Dt / 60); % Dt is in minutes, so divide by 60
    end

    % Store results
    LOL_cumulative_days_all(:, case_num) = LOL_cumulative_hours / 24; % convert hours to days
    LOL_cumulative_minutes_all(:, case_num) = LOL_cumulative_hours * 60; % convert hours to minutes
end

% ---- PLOT 1: Cumulative LOL in Days ----
figure;
plot(time, LOL_cumulative_days_all, 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)', 'FontSize', 16);
ylabel('Cumulative Loss of Life (days)', 'FontSize', 16);
%title('Cumulative Loss of Life in Days');
legend(arrayfun(@(k) sprintf('K = %.1f', k), K_values, 'UniformOutput', false), ...
    'Location', 'northwest', 'FontSize', 12);
ylim([0 2]);
yticks(0:0.5:2)
set(gca, 'FontSize', 14);

% ---- PLOT 2: Cumulative LOL in Minutes ----
figure;
plot(time, LOL_cumulative_minutes_all, 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)');
ylabel('Cumulative Loss of Life (minutes)');
%title('Cumulative Loss of Life in Minutes');
legend(arrayfun(@(k) sprintf('K = %.1f', k), K_values, 'UniformOutput', false), ...
    'Location', 'northwest');
ylim([0 3000]);
yticks(0:500:3000);