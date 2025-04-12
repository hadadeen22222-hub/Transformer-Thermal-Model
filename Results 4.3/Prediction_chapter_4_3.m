clc; clear; close all;

%% -----------------------------
%        Shared Constants
%% -----------------------------
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
Dt = 3;  % Time step in minutes

% --- Correct Initial conditions from 168 hours ---
theta_0_init = 38.27;
delta_theta_h1_init = 31.64;
delta_theta_h2_init = 9.32;

%% -----------------------------
% Setup
%% -----------------------------
T_total = 120;
time = (0:1:T_total-1)';
K_values = 0:0.05:1.4;
K_range = 1:0.01:2;

ambient_temps = [19.2, 22.2, 25.2];
colors = {'b', [0 0.5 1], [0.2 0.7 1]};
line_styles = {'-', '--', ':'};
red_colors = {'r', [1 0.5 0], [1 0.2 0.2]};

%% -----------------------------
% Plot Setup
figure;
yyaxis left
hold on;
legend_entries = {};

% --- Left Y-axis: Final HST ---
for a = 1:length(ambient_temps)
    Tambient = ambient_temps(a);
    final_HST = NaN(1, length(K_values));
    
    for case_num = 1:length(K_values)
        K = ones(length(time), 1) * K_values(case_num);

        TOT = NaN(length(time), 1);
        HST = NaN(length(time), 1);

        TOT(1) = theta_0_init;
        HST(1) = theta_0_init + (delta_theta_h1_init - delta_theta_h2_init);
        delta_theta_h1 = delta_theta_h1_init;
        delta_theta_h2 = delta_theta_h2_init;

        for i = 2:length(time)
            D_theta_0 = (Dt / (k11 * tao_0)) * (((1 + K(i)^2 * R)/(1 + R))^x * delta_theta_or - (TOT(i-1) - Tambient));
            TOT(i) = TOT(i-1) + D_theta_0;

            D_delta_theta_h1 = (Dt / (k22 * tao_w)) * (k21 * delta_theta_hr * K(i)^y - delta_theta_h1);
            delta_theta_h1 = delta_theta_h1 + D_delta_theta_h1;

            D_delta_theta_h2 = (Dt / ((1 / k22) * tao_0)) * ((k21 - 1) * delta_theta_hr * K(i)^y - delta_theta_h2);
            delta_theta_h2 = delta_theta_h2 + D_delta_theta_h2;

            delta_theta_h = delta_theta_h1 - delta_theta_h2;
            HST(i) = TOT(i) + delta_theta_h;
        end

        final_HST(case_num) = HST(end);
    end

    plot(K_values, final_HST, line_styles{a}, 'Color', colors{a}, 'LineWidth', 2);
    legend_entries{end+1} = sprintf('Equilibrium HST (%d°C Ambient)', Tambient);
end

% === Fix left y-axis ===
ylabel('Steady State Hot-Spot Temperature (°C)');
ylim([20 160]);           % Left axis from 20 to 160
yticks(20:20:160);        % Ticks every 20

%% -----------------------------
% Right Y-axis: Time to 120°C
yyaxis right
hold on;

for a = 1:length(ambient_temps)
    Tambient = ambient_temps(a);
    HST_times = NaN(length(K_range), 1);

    for j = 1:length(K_range)
        K = K_range(j);
        theta_0 = theta_0_init;
        delta_theta_h1 = delta_theta_h1_init;
        delta_theta_h2 = delta_theta_h2_init;
        theta_0_prev = theta_0;
        prev_HST = NaN;
        time_elapsed = 0;
        HST_threshold = 120;

        while time_elapsed < 5000
            D_theta_0 = (Dt / (k11 * tao_0)) * (((1 + K^2 * R)/(1 + R))^x * delta_theta_or - (theta_0_prev - Tambient));
            theta_0 = theta_0_prev + D_theta_0;

            D_delta_theta_h1 = (Dt / (k22 * tao_w)) * (k21 * delta_theta_hr * K^y - delta_theta_h1);
            delta_theta_h1 = delta_theta_h1 + D_delta_theta_h1;

            D_delta_theta_h2 = (Dt / ((1 / k22) * tao_0)) * ((k21 - 1) * delta_theta_hr * K^y - delta_theta_h2);
            delta_theta_h2 = delta_theta_h2 + D_delta_theta_h2;

            delta_theta_h = delta_theta_h1 - delta_theta_h2;
            HST = theta_0 + delta_theta_h;

            if HST >= HST_threshold
                if ~isnan(prev_HST)
                    fraction = (HST_threshold - prev_HST) / (HST - prev_HST);
                    HST_times(j) = time_elapsed + fraction * Dt;
                else
                    HST_times(j) = time_elapsed;
                end
                break;
            end

            prev_HST = HST;
            theta_0_prev = theta_0;
            time_elapsed = time_elapsed + Dt;
        end
    end

    plot(K_range, HST_times, line_styles{a}, 'Color', red_colors{a}, 'LineWidth', 2);
    legend_entries{end+1} = sprintf('Time to 120°C (%d°C Ambient)', Tambient);
end

% === Fix right y-axis ===
ylabel('Time to Reach 120°C (minutes)');
ylim([0 175]);            % Right axis from 0 to 175
yticks(0:25:175);         % Ticks every 25

xlabel('Load Factor (K)');

% === Legend and style ===
legend('Steady State HST at 19.2°C','Steady State HST at 22.2°C','Steady State HST at 25.2°C',...
       'Time to safety limit at 19.2°C','Time to safety limit at 22.2°C', 'Time to safety limit at 25.2°C', ...
       'Location', 'northwest');

% === Make tick numbers and labels bigger ===
set(gca, 'FontSize', 15); % Bigger numbers (tick labels)

xlabel('Load Factor (K)', 'FontSize', 17, 'FontWeight', 'bold');

yyaxis left
ylabel('Steady State Hot-Spot Temperature (°C)', 'FontSize', 17, 'FontWeight', 'bold');

yyaxis right
ylabel('Time to Reach 120°C (minutes)', 'FontSize', 17, 'FontWeight', 'bold');

% === Turn grid OFF ===
grid off;