clc; clear; close all;

% Constants
delta_theta_or = 45;     % Top-oil temperature rise at rated losses (°C)
tao_0 = 150;             % Oil time constant (minutes)
R = 8;                   % Ratio of load losses to no-load losses
x = 0.8;                 % Oil exponent
k11 = 0.5;               % Thermal model constant
Tambient = 20;           % Ambient temperature (°C)
Dt = 3;                  % Time step (minutes)

% Time setup
T_total = 120;                          % Total simulation time (minutes)
time = (0:Dt:T_total)';                 % Time vector (aligned with Dt)
time_steps = length(time);             % Number of simulation points

% Load factors to simulate
K_values_TOT = 0.1:0.2:1.7;             % Load factor range
num_cases_TOT = length(K_values_TOT);  % Number of simulations

% Storage for results
TOT_cases = NaN(time_steps, num_cases_TOT);  % Top-oil temperature for each K

% Main simulation loop
for case_num = 1:num_cases_TOT
    K = ones(time_steps, 1) * K_values_TOT(case_num);  % Constant load for each case
    TOT = NaN(time_steps, 1);
    
    % Initial condition (cold start)
    TOT(1) = 20;  % Start at ambient temperature
    
    for i = 2:time_steps
        D_theta_0 = (Dt / (k11 * tao_0)) * ...
            (((1 + K(i)^2 * R) / (1 + R))^x * delta_theta_or - (TOT(i-1) - Tambient));
        TOT(i) = TOT(i-1) + D_theta_0;
    end

    % Store in result matrix
    TOT_cases(:, case_num) = TOT;
end

% Plot Top-Oil Temperature for all K values
figure;
plot(time, TOT_cases, 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)');
ylabel('Top-Oil Temperature (°C)');
legend(arrayfun(@(k) sprintf('K = %.1f', k), K_values_TOT, 'UniformOutput', false), 'Location', 'northwest');
hold off;