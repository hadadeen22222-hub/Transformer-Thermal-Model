clc; clear; close all;

% Constants
delta_theta_or = 45;     % Top-oil temperature rise at rated losses (°C)
tao_0 = 150;             % Oil time constant (minutes)
R = 8;                   % Ratio of load losses to no-load losses
x = 0.8;                 % Oil exponent
k11 = 0.5;               % Thermal model constant
Dt = 3;                  % Time step (minutes)

% Time setup
T_total = 120;                          % Total simulation time (minutes)
time = (0:Dt:T_total)';                 % Time vector (aligned with Dt)
time_steps = length(time);              % Number of simulation points

% Ambient temperatures to simulate
Tambient_values = -10:5:35;             % Ambient temperature range
num_cases_TOT = length(Tambient_values);  % Number of simulations

% Storage for results
TOT_cases = NaN(time_steps, num_cases_TOT);  % Top-oil temperature for each ambient temperature

% Load factor is constant
K = ones(time_steps, 1) * 1;  % Constant load factor K = 1

% Main simulation loop
for case_num = 1:num_cases_TOT
    Tambient = Tambient_values(case_num);  % Current ambient temperature
    TOT = NaN(time_steps, 1);
    
    % Initial condition (cold start at Tambient)
    TOT(1) = Tambient;
    
    for i = 2:time_steps
        D_theta_0 = (Dt / (k11 * tao_0)) * ...
            (((1 + K(i)^2 * R) / (1 + R))^x * delta_theta_or - (TOT(i-1) - Tambient));
        TOT(i) = TOT(i-1) + D_theta_0;
    end

    % Store in result matrix
    TOT_cases(:, case_num) = TOT;
end

% Plot Top-Oil Temperature for all ambient temperatures
figure;
plot(time, TOT_cases, 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)', 'FontSize', 16);
ylabel('Top-Oil Temperature (°C)', 'FontSize', 16);
legend(arrayfun(@(T) sprintf('Tambient = %d°C', T), Tambient_values, 'UniformOutput', false), 'Location', 'northwest','FontSize', 11);
set(gca, 'FontSize', 14);
hold off;


%title('Top-Oil Temperature for Different Ambient Temperatures');