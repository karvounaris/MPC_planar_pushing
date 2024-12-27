%========================================================================%
% Τhis script runs a simulation using the mpc controller
%========================================================================%

%% Set up simulation parameters
clear
close all
clc

% Object's parameters
len = 0.1;
radius = 0.05;
width = radius*2;
height = 0.05;
mass = 3;
rectangular_prism_mass = mass * (2*len*radius) / (2*len*radius + pi*radius^2);
cylinder_mass = mass * (pi*radius^2) / (2*len*radius + pi*radius^2);
object_shape = "rectangular_capsule_prism";
I_object = calculate_inertia_matrix(len, width, height, radius, ...
                                    rectangular_prism_mass, cylinder_mass/2, ...
                                    object_shape);
contact_area = calculate_contact_area(len, width, radius, object_shape);

% Limit surface model
alpha = 0.63;
g = 9.81;
R = sqrt(contact_area/pi);
F_N = mass * g;
mu_ground = 0.5;
mu = 0.5;

L = [1/(mu_ground*F_N)^2 0 0;
     0 1/(mu_ground*F_N)^2 0;
     0 0 1/(alpha*R*mu_ground*F_N)^2];

trajectory_radius = 0.5;

duration = 0.5;
timestep = 0.001;

x_0 = 0;
x_f = 0.06 * duration;
% x_f = 2;
y_0 = 0;
y_f = 0.06 * duration;
% y_f = 3;
% [x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         seventh_trajectory_straight_line(duration, x_0, x_f, y_0, y_f, timestep);

[x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, ~] = ...
                        constant_velocity_trajectory_straight_line(duration, x_0, x_f, y_0, y_f, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         quarter_circle_trajectory(duration, trajectory_radius, timestep);
% x(1,1) = trajectory_radius;

[fn_star, ft_star, phi_star_dot, phi_star, ~] = ...
                          calculate_u_star_one_contact_point(L, x_star_dot, y_star_dot, theta_star, ...
                          theta_star_dot, len, radius, timestep, duration);

x_star = [x_star; y_star; theta_star; phi_star];
u_star = [fn_star; ft_star; phi_star_dot];

% Define parameter grids
Q_values = {
    50 * diag([4, 4, 0.1, 0]);
    80 * diag([4, 4, 0.1, 0]);
    120 * diag([4, 4, 0.1, 0]);
    200 * diag([4, 4, 0.1, 0]);
};
QN_values = {
    15000 * diag([4, 4, 0.1, 0]);
    20000 * diag([4, 4, 0.1, 0]);
    25000 * diag([4, 4, 0.1, 0]);
    30000 * diag([4, 4, 0.1, 0]);
};
R_values = {
    0.005 * diag([1, 1, 0]);
    0.01 * diag([1, 1, 0]);
    0.1 * diag([1, 1, 0])
};
mpc_timesteps = [0.02, 0.03, 0.04, 0.05];

% N_values = [20, 30, 50, 60];
N = 20;

x_initial = [0.05, 0, 0, phi_star(1)];

results = grid_search(x_star, u_star, Q_values, QN_values, R_values, mpc_timesteps,... 
                      N_values, duration, timestep, mu, L, radius, len, mass, x_initial, I_object);

save('./variables.mat');

%% Visualize top 5 results
results_table = struct2table(results);
% % Ensure RMSE is numeric
% if iscell(results_table.RMSE)
%     results_table.RMSE = cell2mat(results_table.RMSE);
% end
% 
% % Sort by RMSE (ascending order)
% sorted_results_table = sortrows(results_table, 'RMSE');
% 
% % Convert Q, QN, and R to string for better visualization
% sorted_results_table.Q_str = cellfun(@(x) mat2str(x, 4), sorted_results_table.Q, 'UniformOutput', false);
% sorted_results_table.QN_str = cellfun(@(x) mat2str(x, 4), sorted_results_table.QN, 'UniformOutput', false);
% sorted_results_table.R_str = cellfun(@(x) mat2str(x, 4), sorted_results_table.R, 'UniformOutput', false);
% 
% % Display the top 5 configurations
% disp(sorted_results_table(1:5, {'Q_str', 'QN_str', 'R_str', 'RMSE', 'TerminalStateDeviation', 'AvgSolverTime'}));

%% Visualize using scatterplot
% Scatter plot of RMSE vs Terminal State Deviation
figure;
scatter(results_table.RMSE, results_table.TerminalStateDeviation, 'filled');
xlabel('RMSE');
ylabel('Terminal State Deviation');
title('RMSE vs Terminal State Deviation');

saveas(gcf,'./figures/Visualize_using_scatterplot.png');

%% Line Plot: RMSE vs. Prediction Horizon
figure;
hold on;
for timestep = mpc_timesteps
    subset = results_table(results_table.mpc_timestep == timestep, :);
    plot(subset.N, subset.RMSE, '-o', 'DisplayName', sprintf('Timestep = %.2f', timestep));
end
xlabel('Prediction Horizon (N)');
ylabel('RMSE');
title('RMSE vs Prediction Horizon');
legend('Location', 'best');
hold off;
saveas(gcf,'./figures/RMSE_vs_Prediction_Horizon.png');

%% Bar Plot: Average Solver Time
% Group data by N and calculate the mean AvgSolverTime
grouped_data = groupsummary(results_table, 'N', 'mean', 'AvgSolverTime');

% Bar plot of AvgSolverTime vs Prediction Horizon
figure;
bar(categorical(grouped_data.N), grouped_data.mean_AvgSolverTime);
xlabel('Prediction Horizon (N)');
ylabel('Average Solver Time (s)');
title('Solver Time vs Prediction Horizon');
saveas(gcf,'./figures/Average_Solver_Time.png');

%% 3D Scatter Plot
figure;
scatter3(results_table.RMSE, results_table.TerminalStateDeviation, results_table.AvgSolverTime, ...
    'filled');
xlabel('RMSE');
ylabel('Terminal State Deviation');
zlabel('Solver Time');
title('3D Scatter Plot of Results');
grid on;
saveas(gcf,'./fiugres/3D_Scatter_Plot.png');

%% RMSE vs Q_values (mean RMSE for each Q_value)
figure;
mean_rmse_Q = [];
for i = 1:length(Q_values)
    % Compare Q values using isequal
    subset = results_table(cellfun(@(x) isequal(x, Q_values{i}), results_table.Q), :);
    mean_rmse_Q = [mean_rmse_Q; mean(subset.RMSE)];
end
bar(1:length(Q_values), mean_rmse_Q);
xticks(1:length(Q_values));
xticklabels(cellfun(@(x) mat2str(x, 4), Q_values, 'UniformOutput', false));
xlabel('Q Values');
ylabel('Mean RMSE');
title('RMSE vs Q Values');
grid on;
saveas(gcf,'./figures/3D_Scatter_Plot.png');

%% RMSE vs QN_values (mean RMSE for each QN_value)
figure;
mean_rmse_QN = [];
for i = 1:length(QN_values)
    % Compare QN values using isequal
    subset = results_table(cellfun(@(x) isequal(x, QN_values{i}), results_table.QN), :);
    mean_rmse_QN = [mean_rmse_QN; mean(subset.RMSE)];
end
bar(1:length(QN_values), mean_rmse_QN);
xticks(1:length(QN_values));
xticklabels(cellfun(@(x) mat2str(x, 4), QN_values, 'UniformOutput', false));
xlabel('QN Values');
ylabel('Mean RMSE');
title('RMSE vs QN Values');
grid on;
saveas(gcf,'./figures/RMSE_vs_QN_values.png');

%% RMSE vs R_values (mean RMSE for each R_value)
figure;
mean_rmse_R = [];
for i = 1:length(R_values)
    % Compare R values using isequal
    subset = results_table(cellfun(@(x) isequal(x, R_values{i}), results_table.R), :);
    mean_rmse_R = [mean_rmse_R; mean(subset.RMSE)];
end
bar(1:length(R_values), mean_rmse_R);
xticks(1:length(R_values));
xticklabels(cellfun(@(x) mat2str(x, 4), R_values, 'UniformOutput', false));
xlabel('R Values');
ylabel('Mean RMSE');
title('RMSE vs R Values');
grid on;
saveas(gcf,'./figures/RMSE_vs_R_values.png');

%% RMSE vs mpc_timesteps (mean RMSE for each timestep)
figure;
mean_rmse_timestep = [];
for i = 1:length(mpc_timesteps)
    % Compare timestep values using direct comparison
    subset = results_table(results_table.mpc_timestep == mpc_timesteps(i), :);
    mean_rmse_timestep = [mean_rmse_timestep; mean(subset.RMSE)];
end
bar(1:length(mpc_timesteps), mean_rmse_timestep);
xticks(1:length(mpc_timesteps));
xticklabels(arrayfun(@(x) sprintf('%.2f', x), mpc_timesteps, 'UniformOutput', false));
xlabel('MPC Timesteps');
ylabel('Mean RMSE');
title('RMSE vs MPC Timesteps');
grid on;
saveas(gcf,'./figures/RMSE_vs_mpc_timesteps.png');

%% RMSE vs N_values (mean RMSE for each N_value)
figure;
mean_rmse_N = [];
for i = 1:length(N_values)
    % Compare N values using direct comparison
    subset = results_table(results_table.N == N_values(i), :);
    mean_rmse_N = [mean_rmse_N; mean(subset.RMSE)];
end
bar(1:length(N_values), mean_rmse_N);
xticks(1:length(N_values));
xticklabels(arrayfun(@(x) sprintf('%d', x), N_values, 'UniformOutput', false));
xlabel('Prediction Horizon (N)');
ylabel('Mean RMSE');
title('RMSE vs Prediction Horizon (N)');
grid on;
saveas(gcf,'./figures/RMSE_vs_N_values.png');
