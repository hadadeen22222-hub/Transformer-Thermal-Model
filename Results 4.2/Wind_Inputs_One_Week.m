% Load the wind speed data
wind_speed_data = readmatrix('Wind Speed.csv');
wind_inputs_data = readmatrix('Wind_Inputs11.csv');

% Define original time step (10 minutes)
Dt_original = 10; % 10-minute intervals
time_original = (0:Dt_original:(length(wind_speed_data) - 1) * Dt_original)';

% Define new time step (1 minute)
Dt_new = 1; % 1-minute intervals
time_new = (0:Dt_new:(time_original(end)))';

% Perform linear interpolation
wind_speed_interpolated = interp1(time_original, wind_speed_data, time_new, 'linear');
wind_inputs_interpolated = interp1(time_original, wind_inputs_data, time_new, 'linear');

% Convert time to days for plotting
time_days_new = time_new / (24 * 60);

% Define the duration for one week (in minutes)
one_week_minutes = 7 * 24 * 60;

% Extract the first week's data for wind speed
indices_one_week = find(time_new <= one_week_minutes);
time_days_week = time_days_new(indices_one_week);
wind_speed_week = wind_speed_interpolated(indices_one_week);
wind_inputs_week = wind_inputs_interpolated(indices_one_week);

% Reduce density by selecting every 50th data point
downsample_factor_week = 50;
indices_reduced_week = 1:downsample_factor_week:length(time_days_week);

% Plot wind speed for one week
figure('Position', [100, 100, 1200, 600]);
plot(time_days_week(indices_reduced_week), wind_speed_week(indices_reduced_week), '-b', 'LineWidth', 1.5);

grid on;
xlabel('Time (Days)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Wind Speed (m/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Wind Speed Over One Week', 'FontSize', 14, 'FontWeight', 'bold');
legend('Wind Speed (m/s)', 'Location', 'best', 'FontSize', 12);

% Filter out negative input values
wind_inputs_week(wind_inputs_week < 0) = 0;

% Plot wind farm inputs for one week
figure('Position', [100, 100, 1200, 600]);
plot(time_days_week(indices_reduced_week), wind_inputs_week(indices_reduced_week), '-r', 'LineWidth', 1.5);

grid on;
xlabel('Time (Days)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Input Value (Per Unit)', 'FontSize', 18, 'FontWeight', 'bold');
legend('Wind Farm Inputs', 'Location', 'best', 'FontSize', 16, 'FontWeight', 'bold');

ylim([0 inf]);  % Set y-axis to start at 0

% Increase tick label font size
set(gca, 'FontSize', 20);