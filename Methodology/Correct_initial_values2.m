% Semester 1 model BUT with correct initial conditions equations 

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

% Input Data
K = [0.81, 0.87, 0.88, 0.86, 0.90, 0.92, 0.95, 0.96, 0.97, 1.00, ...
     1.70, 1.70, 1.73, 1.72, 1.69, 1.68, 1.71, 1.69, 1.67, 1.68, ...
     1.63, 1.59, 1.53, 1.49, 1.41, 1.38, 1.32, 1.28, 1.21, 1.19, ...
     0.87, 0.88, 0.87, 0.86, 0.85, 0.87, 0.83, 0.86, 0.85, 0.82, ...
     0.86];

theta_a = [30.3, 29.9, 29.8, 29.5, 29.6, 29.5, 29.5, 28.9, 29.0, 28.6, ...
           28.0, 28.7, 27.8, 28.1, 27.9, 27.1, 26.9, 26.7, 27.2, 26.7, ...
           26.9, 26.5, 26.2, 26.3, 25.4, 25.6, 25.3, 24.8, 24.5, 24.3, ...
           24.1, 24.3, 24.1, 23.4, 23.6, 23.8, 23.1, 23.3, 23.1, 22.3, ...
           22.2];

% Time step
Dt = 3;  % Minutes
time_steps = length(K);
time = (0:Dt:(time_steps - 1) * Dt)'; % Time array

% Initialize arrays
HST = NaN(time_steps, 1); % Hot-spot temperature
TOT = NaN(time_steps, 1); % Top-oil temperature
PU_life = NaN(time_steps, 1); % Per-unit life
LOL = 0; % Total loss of life

% Initial conditions
K_0 = K(1);
theta_a_0 = theta_a(1);
theta_0 = ((1+K_0.^2.*R)./(1+R)).^x .* delta_theta_or + theta_a_0;  % Initial top-oil temperature
delta_theta_h1 = k21 * K_0.^y * delta_theta_hr;
delta_theta_h2 = (k21 - 1) * K_0.^y * delta_theta_hr;

% Solve iteratively
for i = 1:time_steps
    % Top-oil temperature
    D_theta_0 = (Dt / (k11 * tao_0)) * (((1 + K(i)^2 * R) / (1 + R))^x * delta_theta_or - (theta_0 - theta_a(i)));
    theta_0 = theta_0 + D_theta_0;

    % Winding gradients
    D_delta_theta_h1 = (Dt / (k22 * tao_w)) * (k21 * delta_theta_hr * K(i)^y - delta_theta_h1);
    delta_theta_h1 = delta_theta_h1 + D_delta_theta_h1;

    D_delta_theta_h2 = (Dt / ((1 / k22) * tao_0)) * ((k21 - 1) * delta_theta_hr * K(i)^y - delta_theta_h2);
    delta_theta_h2 = delta_theta_h2 + D_delta_theta_h2;

    % Hot-spot gradient
    delta_theta_h = delta_theta_h1 - delta_theta_h2;

    % Store results
    TOT(i) = theta_0;
    HST(i) = theta_0 + delta_theta_h;

    % Per-unit life
    PU_life(i) = exp(15000 / 383 - 15000 / (HST(i) + 273));

    % Accumulate loss of life
    LOL = LOL + (PU_life(i) * Dt / 60); % Convert to hours

    % Cumulative LOL
    LOL_cumulative_minutes(i) = sum(PU_life(1:i) * Dt);

end

% Convert cumulative LOL from minutes to days
LOL_cumulative_days = LOL_cumulative_minutes / (24 * 60);

% Create table for results
results_table = table(time, TOT, HST, PU_life, 'VariableNames', ...
    {'Time_Minutes', 'Top_Oil_Temperature_C', 'Hot_Spot_Temperature_C', 'Per_Unit_Life'});

disp(results_table);

%% Plots

figure;
plot(time, HST,'-r', 'LineWidth', 1.5); hold on;
grid on;
xlabel('Time (minutes)', 'FontSize', 15);
ylabel('Temperature (°C)', 'FontSize', 15);
legend('Hot-Spot Temperature', 'FontSize', 12);
set(gca, 'FontSize', 14);

figure;
plot(time, TOT, '-b', 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)', 'FontSize', 15);
ylabel('Temperature (°C)', 'FontSize', 15);
legend('Top-Oil Temperature', 'FontSize', 12);
set(gca, 'FontSize', 14);

figure;
plot(time, LOL_cumulative_minutes, '-m', 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)', 'FontSize', 15);
ylabel('Cumulative Loss of Life (minutes)', 'FontSize', 15);
set(gca, 'FontSize', 14);

figure;
plot(time, LOL_cumulative_days, '-g', 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)', 'FontSize', 15);
ylabel('Cumulative Loss of Life (days)', 'FontSize', 15);
legend('Cumulative Loss of Life', 'FontSize', 12);
set(gca, 'FontSize', 14);

figure;
plot(time, HST, '-rx', 'LineWidth', 1.5); hold on;
plot(time, HST, '-b', 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)', 'FontSize', 15);
ylabel('Temperature (°C)', 'FontSize', 15);
legend('HST (IEC)', 'HST (MATLAB)', 'FontSize', 14);
set(gca, 'FontSize', 14);

figure;
plot(time, TOT, '-rx', 'LineWidth', 1.5); hold on;
plot(time, TOT, '-b', 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)', 'FontSize', 15);
ylabel('Temperature (°C)', 'FontSize', 15);
legend('TOT (IEC)', 'TOT (MATLAB)', 'FontSize', 14);
set(gca, 'FontSize', 14);

figure;
plot(time, LOL_cumulative_days, '-rx', 'LineWidth', 1.5); hold on;
plot(time, LOL_cumulative_days, '-b', 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)', 'FontSize', 15);
ylabel('Cumulative Loss of Life (days)', 'FontSize', 15);
legend('LOL (IEC)', 'LOL (MATLAB)', 'FontSize', 14);
set(gca, 'FontSize', 14);


figure;
plot(time, K, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)', 'FontSize', 16);
ylabel('K Values', 'FontSize', 16);
set(gca, 'FontSize', 14);

figure;
plot(time, theta_a, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)', 'FontSize', 16);
ylabel('Ambient Temperature (°C)', 'FontSize', 16);
set(gca, 'FontSize', 14);

figure;
yyaxis left
plot(time, theta_a, '-o', 'LineWidth', 1.5);
ylabel('Ambient Temperature (°C)', 'FontSize', 16);

yyaxis right
plot(time, K, '-o', 'LineWidth', 1.5);
ylabel('K Values', 'FontSize', 16);

xlabel('Time (minutes)', 'FontSize', 16);
title('Ambient Temperature and Load Factor Over Time', 'FontSize', 18);
legend({'Ambient Temperature', 'Load Factor'}, 'Location', 'northwest', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 14);

% Merged plot for HST and TOT
figure;
plot(time, HST, '-r', 'LineWidth', 1.5); hold on;
plot(time, TOT, '-b', 'LineWidth', 1.5); 
grid on;
xlabel('Time (minutes)', 'FontSize', 15);
ylabel('Temperature (°C)', 'FontSize', 15);
legend({'Hot-Spot Temperature', 'Top-Oil Temperature'}, 'FontSize', 12, 'Location', 'best');
%title('Hot-Spot and Top-Oil Temperatures Over Time', 'FontSize', 16);
set(gca, 'FontSize', 14);

% Single merged plot for HST and TOT (with IEC and MATLAB markers)
figure;
plot(time, HST, '-rx', 'LineWidth', 1.5); hold on;    % HST IEC
plot(time, HST, 'b', 'LineWidth', 1.5);              % HST MATLAB
plot(time, TOT, '-rx', 'LineWidth', 1.5);             % TOT IEC
plot(time, TOT, 'b', 'LineWidth', 1.5);              % TOT MATLAB
grid on;
xlabel('Time (minutes)', 'FontSize', 15);
ylabel('Temperature (°C)', 'FontSize', 15);
legend({'HST (IEC)', 'HST (MATLAB)', 'TOT (IEC)', 'TOT (MATLAB)'}, 'FontSize', 12, 'Location', 'best');
%title('Hot-Spot and Top-Oil Temperature Comparison', 'FontSize', 16);
set(gca, 'FontSize', 14);