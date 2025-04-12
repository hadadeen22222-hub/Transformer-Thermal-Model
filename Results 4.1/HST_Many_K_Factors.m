clc; clear; close all;

% Constants
delta_theta_or = 45;     % Top-oil temperature rise at rated losses, K
delta_theta_hr = 35;     % Hot-spot-to-top-oil gradient at rated current, K
tao_0 = 150;             % Average oil time constant, min
tao_w = 7;               % Winding time constant, min
R = 8;                   % Ratio of load losses at rated current to no-load losses
x = 0.8;                 % Oil exponent
y = 1.3;                 % Winding exponent
k11 = 0.5;               % Thermal model constant
k21 = 2;                 % Thermal model constant
k22 = 2;                 % Thermal model constant
Tambient = 20;           % Initial ambient temperature set to 20°C (Cold Start)
Dt = 3;                  % Time step, minutes

% Time range
T_total = 120;  % Simulation time (0 to 120 min)
time = (0:1:T_total-1)'; % Time vector (0-120 minutes)

% Define load factors
K_values = 0.1:0.2:1.7;  % Load factor K from 0.1 to 1.9 in steps of 0.2
num_cases = length(K_values);

% Initialize storage variables
TOT_cases = NaN(length(time), num_cases);
HST_cases = NaN(length(time), num_cases);
PU_life_cases = NaN(length(time), num_cases);
LOL_cumulative_cases = NaN(length(time), num_cases);
final_HST = NaN(1, num_cases); % Store final hot-spot temperature for each K

% Run simulations for each K value
for case_num = 1:num_cases
    K = ones(length(time), 1) * K_values(case_num);
    
    % Initialize temperature and life loss variables
    TOT = NaN(length(time), 1);
    HST = NaN(length(time), 1);
    
    % Initial conditions (Cold Start at 20°C)
    TOT(1) = 20;
    HST(1) = 20;
    
    delta_theta_h1 = 0;
    delta_theta_h2 = 0;
    
    % Iteration loop
    for i = 2:length(time)
        D_theta_0 = (Dt / (k11 * tao_0)) * (((1 + K(i)^2 * R) / (1 + R))^x * delta_theta_or - (TOT(i-1) - Tambient));
        TOT(i) = TOT(i-1) + D_theta_0;

        D_delta_theta_h1 = (Dt / (k22 * tao_w)) * (k21 * delta_theta_hr * K(i)^y - delta_theta_h1);
        delta_theta_h1 = delta_theta_h1 + D_delta_theta_h1;

        D_delta_theta_h2 = (Dt / ((1 / k22) * tao_0)) * ((k21 - 1) * delta_theta_hr * K(i)^y - delta_theta_h2);
        delta_theta_h2 = delta_theta_h2 + D_delta_theta_h2;

        delta_theta_h = delta_theta_h1 - delta_theta_h2;
        HST(i) = TOT(i) + delta_theta_h;
    end
    
    % Store results
    TOT_cases(:, case_num) = TOT;
    HST_cases(:, case_num) = HST;
    final_HST(case_num) = HST(end); % Store the final HST value for this K
end

% Plot Hot-Spot Temperature Comparison
figure;
plot(time, HST_cases, 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)');
ylabel('Temperature (°C)');
legend(arrayfun(@(k) sprintf('K=%.1f', k), K_values, 'UniformOutput', false));
hold off;

% Plot Final Hot-Spot Temperature vs Load Factor
figure;
plot(K_values, final_HST, 'o-', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [1, 0.5, 0]); % Orange color
grid on;
xlabel('Load Factor (K)');
ylabel('Final Hot-Spot Temperature (°C)');
hold off;