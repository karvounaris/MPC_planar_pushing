%========================================================================%
% Τhis script runs a simulation using the mpc controller
%========================================================================%

%% Set up simulation parameters
clear
close all
clc

import casadi.*

% len = 0.2;
% wid = 0.15;
% radius = 0.075;
% height = 0.18;
% mass = 1.1;
% object_shape = "rectangular_prism";
len = 0.1;
radius = 0.05;
wid = radius*2;
height = 0.05;
mass = 4;
object_shape = "rectangular_capsule_prism";
rectangular_prism_mass = mass * (2*len*radius) / (2*len*radius + pi*radius^2);
cylinder_mass = mass * (pi*radius^2) / (2*len*radius + pi*radius^2);
I_object = calculate_inertia_matrix(len, wid, height, radius, ...
                                    rectangular_prism_mass, cylinder_mass/2, mass, ...
                                    object_shape);
contact_area = calculate_contact_area(len, wid, radius, object_shape);

% Limit surface model
alpha = 0.63;
g = 9.81;
R = sqrt(contact_area/pi);
F_N = mass * g;
mu_ground = 0.5;
mu = 0.2;

L = [1/(mu_ground*F_N)^2 0 0;
     0 1/(mu_ground*F_N)^2 0;
     0 0 1/(alpha*R*mu_ground*F_N)^2];

% Izz_rectangle = (1/12)*len*wid^3 + (1/12)*wid*len^3;
% Izz_halfdisks = (pi/2)*radius^4;
% A_halfdisk = (pi/2)*radius^2;
% Izz_shift = A_halfdisk * (len/2)^2 * 2;
% Izz_patch = Izz_rectangle + Izz_halfdisks + Izz_shift;
% I_object(3,3) = Izz_patch;
% 
% f_max = mu_ground * mass * g;
% m_max = f_max * (Izz_patch / contact_area);
% 
% L = diag([1/f_max^2, 1/f_max^2, 1/m_max^2]);

x_0 = 0;
y_0 = 0;
x_f = 0.3;
y_f = 0.3;
v_constant = 0.05;
timestep = 0.002;
trajectory_radius = 0.2;

v_constant_s = 0.05;
x_center = x_0 - trajectory_radius;
y_center = y_0;

mpc_timestep = 0.04;
N = 30;

[x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
                        constant_velocity_trajectory_straight_line(v_constant, x_0, x_f, y_0, y_f, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                     constant_velocity_quarter_circle_trajectory(trajectory_radius, v_constant_s, timestep, x_center, y_center);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                     constant_velocity_semi_circle_trajectory(trajectory_radius, v_constant_s, timestep, x_center, y_center);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                     constant_velocity_s_shape_trajectory(trajectory_radius, v_constant_s, timestep, x_center, y_center);

[fn_star, ft_star, phi_star_dot, phi_star] = ...
                          calculate_u_star_one_contact_point(L, x_star_dot, y_star_dot, theta_star, ...
                          theta_star_dot, len, radius, timestep, duration, wid, object_shape);

x_star = [x_star; y_star; theta_star; phi_star];
u_star = [fn_star; ft_star; phi_star_dot];
for i = 1:length(u_star(3,:))
    if u_star(3,i) ~= 0
        u_star(3,i) = 0;
    end
end
% Define the extension for x_star and u_star
x_star_extension = repmat(x_star(:, end), 1, mpc_timestep/timestep + N*mpc_timestep/timestep);
u_star_extension = zeros(size(u_star, 1), mpc_timestep/timestep + N*mpc_timestep/timestep);
x_star = [x_star, x_star_extension];
u_star = [u_star, u_star_extension];

% System's parameters initialization
x = [0; 0; 0; 0];
x_dot = [0; 0; 0; 0];
x_ddot = [0; 0; 0];
u = [0; 0; 0];
x(:,1) = [x_0 + 0.05, y_0 - 0.05, 0, 3*pi/2];
% obstacle_data = [[0.1; 0.2; 0.05], [0.5; 0.1; 0.05], [0.3; 0.4; 0.05]];  % [x_obstacle, y_obstacle, radius_obstacle]
% obstacle_data = [-0.12; 0.04; 0.04];  % [x_obstacle, y_obstacle, radius_obstacle]
obstacle_data = [];
n_obstacles = size(obstacle_data, 2);

% % MPC controller tunable parameters
% Q = 80 * diag([10, 10, 0.1, 0]);
% QN = 42000 * diag([10, 10, 0.1, 0]);
% R = 0.02 * diag([1, 1, 0.01]);
% W = 0.01 * diag([0.95, 1, 1]);

% MPC controller tunable parameters
Q = 80 * diag([10, 10, 0.01, 0]);
QN = 42000 * diag([10, 10, 0.1, 0]);
R = 0.01 * diag([1, 1, 0.01]);
W = 0.01 * diag([0.95, 1, 1]);

%% Run simulation
% profile on
dp = [0; 0; 0];
time(1) = 0;
mpc_timestamps = 0;
solver_times = 0;
gurobi_solve_times = 0;
last_time_MIQP_start = 0;

X_init = repmat(x(:,1), 1, N+1);
U_init = ones(3, N);
Z_init = zeros(3, 4);

[solver_controller, args] = create_MPC_MIQP_controller( ...
    Q, QN, R, W, N, mpc_timestep, mu, L, radius, len, wid, object_shape, n_obstacles);

args.x0 = [reshape(X_init, 4*(N+1),1); reshape(U_init, 3*N,1); reshape(Z_init, 3*4,1)];

ground_friction = zeros(3,1);
is_start = 1;
k = 1;

for i = 1:floor(duration/timestep)

    dx(:,i) = x(:,i) - x_star(:,i);
    
    if i == 1
        last_time_MIQP_start = i * timestep;

        for j = 1:N+1
            x_star_mpc(:,j) = x_star(:, i + round((j-1)*mpc_timestep/timestep));
            u_star_mpc(:,j) = u_star(:, i + round((j-1)*mpc_timestep/timestep));
        end

        X0_val = x(:,i);
        X_star_val = x_star_mpc;
        U_star_val = u_star_mpc(:, 1:end-1);
        obstacles_val = obstacle_data;
        
        param_val = [X0_val; X_star_val(:); U_star_val(:); obstacles_val(:)];

        tic;
        sol = solver_controller('p', param_val, ...
                                'x0', args.x0, ...
                                'lbx', args.lbx, 'ubx', args.ubx, ...
                                'lbg', args.lbg, 'ubg', args.ubg);

        solver_times(k) = round(toc*1000) / 1000;
        info = solver_controller.stats();
        gurobi_solve_times(k) = round(info.t_wall_total * 1000) / 1000;
        mpc_timestamps(k) = time(i);
        fprintf('Time is: %g\n', time(i));
        fprintf('Error is: %2.4f %2.4f %2.4f %2.4f %2.4f\n', dx(1,i), dx(2,i), dx(3,i), dx(4,i), ... 
                     sqrt(dx(1,i)^2 + dx(2,i)^2));

        X_opt = reshape(sol.x(1:4*(N+1)), 4, N+1);
        U_opt = reshape(sol.x(4*(N+1)+1:4*(N+1)+3*N), 3, N);
        Z_opt = reshape(sol.x(4*(N+1)+3*N+1:end), 3, 4);

        X_init = [X_opt(:,2:end), X_opt(:,end)];
        U_init = [U_opt(:,2:end), U_opt(:,end)];
        Z_init = [Z_opt(:,2:end), Z_opt(:,end)];

        args.x0 = [reshape(X_init, 4*(N+1),1); reshape(U_init, 3*N,1); reshape(Z_init, 3*4,1)];
        u(:,i) = full(U_opt(:,1));
        z(:,i) = full(Z_opt(:,1));
        k = k + 1;
    elseif (i * timestep - last_time_MIQP_start > solver_times(end) || i == 2)
        last_time_MIQP_start = i * timestep;
        for j = 1:N+1
            x_star_mpc(:,j) = x_star(:, i + round((j-1)*mpc_timestep/timestep));
            u_star_mpc(:,j) = u_star(:, i + round((j-1)*mpc_timestep/timestep));
        end

        u(:,i) = full(U_opt(:,1));
        z(:,i) = full(Z_opt(1)); 

        X0_val = x(:,i);
        X_star_val = x_star_mpc;
        U_star_val = u_star_mpc(:, 1:end-1);
        obstacles_val = obstacle_data;
        
        param_val = [X0_val; X_star_val(:); U_star_val(:); obstacles_val(:)];

        tic;
        sol = solver_controller('p', param_val, ...
                                'x0', args.x0, ...
                                'lbx', args.lbx, 'ubx', args.ubx, ...
                                'lbg', args.lbg, 'ubg', args.ubg);
        
        solver_times(k) = round(toc*1000) / 1000;
        info = solver_controller.stats();
        gurobi_solve_times(k) = round(info.t_wall_total * 1000) / 1000;
        mpc_timestamps(k) = time(i);
        fprintf('Time is: %g\n', time(i));
        fprintf('Error is: %2.4f %2.4f %2.4f %2.4f %2.4f\n', dx(1,i), dx(2,i), dx(3,i), dx(4,i), ... 
                     sqrt(dx(1,i)^2 + dx(2,i)^2));

        if info.unified_return_status == "SOLVER_RET_SUCCESS"
            X_opt = reshape(sol.x(1:4*(N+1)), 4, N+1);
            U_opt = reshape(sol.x(4*(N+1)+1:4*(N+1)+3*N), 3, N);
            Z_opt = reshape(sol.x(4*(N+1)+3*N+1:end), 3, 4);
    
            X_init = [X_opt(:,2:end), X_opt(:,end)];
            U_init = [U_opt(:,2:end), U_opt(:,end)];
            Z_init = [Z_opt(:,2:end), Z_opt(:,end)];
        else
            X_init = [X_opt(:,2:end), X_opt(:,end)];
            U_init = [U_opt(:,2:end), U_opt(:,end)];
            Z_init = [Z_opt(:,2:end), Z_opt(:,end)];
            args.x0 = [reshape(X_init, 4*(N+1),1); reshape(U_init, 3*N,1); reshape(Z_init, 3*4,1)];
        end
        u(:,i) = full(U_opt(:,1));
        z(:,i) = full(Z_opt(:,1));
        k = k + 1;
    elseif i * timestep - last_time_MIQP_start > mpc_timestep
        X_opt = X_init;
        U_opt = U_init;
        Z_opt = Z_init;

        X_init = [X_opt(:,2:end), X_opt(:,end)];
        U_init = [U_opt(:,2:end), U_opt(:,end)];
        Z_init = [Z_opt(:,2:end), Z_opt(:,end)];
        args.x0 = [reshape(X_init, 4*(N+1),1); reshape(U_init, 3*N,1); reshape(Z_init, 3*4,1)];

        u(:,i) = full(U_opt(:,1));
        z(:,i) = full(Z_opt(:,1));
    else
        u(:,i) = u(:,i-1);
        z(:,i) = z(:,i-1);
    end
    
    % Calculate parameters for the motion equation
    w = calculate_motion_model_parameters(u(1:2,i), x(3,i), len, radius, x(4,i), wid, object_shape);
    ground_friction_parameter = 1;
    [gr_frict, number] = calculate_friction_with_ground(L, dp(:,i), ground_friction_parameter);

    if number == 1
        j_fric(i) = 1;
    elseif number == 0
        j_fric(i) = 0;
    end

    ground_friction(:,i) = -gr_frict;
    wrench(:,i) = w;

    x_dot(4,i+1) = u(3,i);
    x(4, i+1) = x(4, i) + x_dot(4, i+1) * timestep;

    % x_ddot(1:3, i+1) = inv(diag([mass mass I_object(3,3)])) * (-gr_frict + w);
    x_ddot(1:3, i+1) = diag([mass mass I_object(3,3)]) \ (-gr_frict + w);
    x_dot(1:3, i+1) = x_dot(1:3, i) + x_ddot(1:3, i+1) .* timestep;
    x(1:3, i+1) = x(1:3, i) + x_dot(1:3, i+1) .* timestep;

    time(i+1) = time(i) + timestep;
    dp(:,i+1) = [x_dot(1, i+1); x_dot(2, i+1); x_dot(3, i+1)];
    
end
% profile viewer



%% Plot Capsule Shape Along Trajectory with Contact Points
figure;
plot(x(1,:), x(2,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Trajectory');
% title('Trajectory with Capsule Shape');
xlabel('x (m)');
ylabel('y (m)');
grid on;
axis equal;

hold on;

plot(x(1,1), x(2,1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(x(1,end), x(2,end), 'yo', 'MarkerFaceColor', 'y', 'DisplayName', 'End');

if object_shape == "rectangular_capsule_prism"
    shape_handle_handle = [];
    % Plot the object shape at several points along the trajectory
    for i = 1:round(length(x(1,:))/10):length(x(1,:))
        shape_handle = get_capsule_shape(len, radius, x(1,i), x(2,i), x(3,i));
    
        if isempty(shape_handle_handle) 
            shape_handle_handle = fill(shape_handle(:,1), shape_handle(:,2), 'b', 'FaceAlpha', 0.2, 'DisplayName', 'Capsule Shape');
        else
            fill(shape_handle(:,1), shape_handle(:,2), 'b', 'FaceAlpha', 0.2);
        end
    end

elseif object_shape == "rectangular_prism"
    shape_handle_handle = [];
    % Plot the object shape at several points along the trajectory
    for i = 1:round(length(x(1,:))/10):length(x(1,:))
        shape_handle = get_rectangle_shape(len, radius, x(1,i), x(2,i), x(3,i));
    
        if isempty(shape_handle_handle) 
            shape_handle_handle = fill(shape_handle(:,1), shape_handle(:,2), 'b', 'FaceAlpha', 0.2, 'DisplayName', 'Capsule Shape');
        else
            fill(shape_handle(:,1), shape_handle(:,2), 'b', 'FaceAlpha', 0.2);
        end
    end
end

contact_x_world = [];
contact_y_world = [];

% Loop to plot contact points and unit vectors separately
for i = 1:length(x(1,:))
    [x_c, y_c, ~, n_c, t_c] = calculate_r_c(x(4,i), len, radius, wid, object_shape);

    R = [cos(x(3,i)), -sin(x(3,i)); sin(x(3,i)), cos(x(3,i))];
    contact_point_global = R * [x_c; y_c] + [x(1,i); x(2,i)];

    contact_x_world = [contact_x_world; contact_point_global(1)];
    contact_y_world = [contact_y_world; contact_point_global(2)];

    if i ~=1
        contact_x_dot_world(i) = (contact_x_world(i) - contact_x_world(i-1))/timestep;
        contact_y_dot_world(i) = (contact_y_world(i) - contact_y_world(i-1))/timestep;
        contact_x_dot_body(i) = cos(x(3)) * contact_x_dot_world(i) - sin(x(3)) * contact_y_dot_world(i);
        contact_y_dot_body(i) = sin(x(3)) * contact_x_dot_world(i) + cos(x(3)) * contact_y_dot_world(i);
    end
end

% Logical indices for j == 1 and j == 0
idx_1 = (j == 1);
idx_0 = (j == 0);

% Plot for j == 1 (red circles)
% plot(contact_x_world(idx_1), contact_y_world(idx_1), 'ro', 'MarkerSize', 2, 'DisplayName', 'Contact Point');
contact_handle = plot(contact_x_world, contact_y_world, 'r-', 'LineWidth', 2, 'DisplayName', 'Contact Point');

% Plot for j == 0 (green circles)
% plot(contact_x_world(idx_0), contact_y_world(idx_0), 'go', 'MarkerSize', 2);
% plot(contact_x_world(idx_0), contact_y_world(idx_0), 'ro', 'MarkerSize', 2);

% plot(contact_x, contact_y, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Contact Path');

% Plot x_star and y_star trajectory
plot(x_star(1,:), x_star(2,:), 'k-', 'LineWidth', 2, 'DisplayName', 'Desired trajectory'); % x_star trajectory in black
plot(x_star(1,1), x_star(2,1), 'go', 'MarkerFaceColor', 'g'); % Start point of x_star in green
plot(x_star(1,end), x_star(2,end), 'yo', 'MarkerFaceColor', 'y'); % End point of x_star in yellow

% Plot the obstacle as a circle
obstacle_handles = []; % store all handles for legend

for j = 1 : size(obstacle_data, 2)
    obs_x = obstacle_data(1, j);
    obs_y = obstacle_data(2, j);
    obs_r = obstacle_data(3, j);

    obs_shape = get_circle_shape(obs_r, obs_x, obs_y);

    if j == 1
        % Add DisplayName only for first (for legend)
        h = fill(obs_shape(:,1), obs_shape(:,2), 'k', ...
            'FaceAlpha', 0.3, 'EdgeColor', 'k', 'LineWidth', 1.5, ...
            'DisplayName', 'Obstacle');
    else
        % No DisplayName for duplicates
        h = fill(obs_shape(:,1), obs_shape(:,2), 'k', ...
            'FaceAlpha', 0.3, 'EdgeColor', 'k', 'LineWidth', 1.5);
    end

    obstacle_handles = [obstacle_handles; h];
end

if n_obstacles >=1
    legend([shape_handle_handle; findobj(gca, 'DisplayName', 'Trajectory'); ...
            findobj(gca, 'DisplayName', 'Start'); findobj(gca, 'DisplayName', 'End'); ...
            findobj(gca, 'DisplayName', 'Desired trajectory');
            contact_handle; obstacle_handles(1)], 'Location', 'Best');
else
    legend([shape_handle_handle; findobj(gca, 'DisplayName', 'Trajectory'); ...
        findobj(gca, 'DisplayName', 'Start'); findobj(gca, 'DisplayName', 'End'); ...
        contact_handle], 'Location', 'Best');
end

hold off;
%% Plot x-x_star, y-y_star, theta-theta_star and phi-phi_star

figure;

subplot(4, 1, 1);
plot(time, x(1, :), 'b-', 'LineWidth', 2);
hold on;
plot(time, x_star(1, 1:length(x(1,:))), 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('x (m)');
legend('x', 'x^*', 'Location', 'best');
grid on;

subplot(4, 1, 2);
plot(time, x(2, :), 'b-', 'LineWidth', 2);
hold on;
plot(time, x_star(2, 1:length(x(2,:))), 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('y (m)');
legend('y', 'y^*', 'Location', 'best');
grid on;

subplot(4, 1, 3);
plot(time, x(3, :), 'b-', 'LineWidth', 2);
hold on;
plot(time, x_star(3, 1:length(x(3,:))), 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('\theta (rad)');
legend('\theta', '\theta^*', 'Location', 'best');
grid on;

subplot(4, 1, 4);
plot(time, x(4, :), 'b-', 'LineWidth', 2);
hold on;
plot(time, x_star(4, 1:length(x(4,:))), 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\phi$ (rad)', 'Interpreter', 'latex');
legend('\phi', '\phi^*', 'Location', 'best');
grid on;

%% Plot x_dot-x_star_dot, y_dot-y_star_dot, theta_dot-theta_star_dot and phi_dot-phi_star_dot

figure;

subplot(4, 1, 1);
plot(time, x_dot(1,:), 'b-', 'LineWidth', 2);
hold on;
plot(time, x_star_dot(1,:), 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{x}$ (m/s)', 'Interpreter', 'latex');
legend({'$\dot{x}$', '$\dot{x}^*$'}, 'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(4, 1, 2);
plot(time, x_dot(2,:), 'b-', 'LineWidth', 2);
hold on;
plot(time, y_star_dot, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{y}$ (m/s)', 'Interpreter', 'latex');
legend({'$\dot{y}$', '$\dot{y}^*$'}, 'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(4, 1, 3);
plot(time, x_dot(3,:), 'b-', 'LineWidth', 2);
hold on;
plot(time, theta_star_dot, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter', 'latex');
legend({'$\dot{\theta}$', '$\dot{\theta}^*$'}, 'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(4, 1, 4);
plot(time, x_dot(4,:), 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{\phi}$ (rad/s)', 'Interpreter', 'latex');
legend({'$\dot{\phi}$'}, 'Interpreter', 'latex', 'Location', 'best');
grid on;

%% Plot dx dy dtheta and dphi in seperate plots

figure;
subplot(4, 1, 1);
plot(time(1:end-1), dx(1,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\bar{x}$ (m)', 'Interpreter', 'latex');
grid on;

subplot(4, 1, 2);
plot(time(1:end-1), dx(2,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\bar{y}$ (m)', 'Interpreter', 'latex');
grid on;

subplot(4, 1, 3);
plot(time(1:end-1), dx(3,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\bar{\theta}$ (rad)', 'Interpreter', 'latex');
grid on;

subplot(4, 1, 4);
plot(time(1:end-1), dx(4,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\bar{\phi}$ (rad)', 'Interpreter', 'latex');
grid on;

%% Plot fn ft, phi_dot_plus and phi_dot_minus in seperate plots

figure;
subplot(3, 1, 1);
plot(time(1:end-1), u(1,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('fn (N)');
grid on;

subplot(3, 1, 2);
plot(time(1:end-1), u(2,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('ft (N)');
grid on;

subplot(3, 1, 3);
plot(time(1:end-1), u(3,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{\phi3}$ (rad/sec)', 'Interpreter', 'latex');
grid on;


%% Plot z1 z2 and z3 in seperate plots

figure;
% First subplot for z1 over time
subplot(3, 1, 1);
plot(time(1:end-1), z(1,:), 'LineWidth', 2);
% title('z1 mode (sticking) Over Time');
xlabel('Time (s)');
ylabel('z1');
grid on;

% Second subplot for z2 over time
subplot(3, 1, 2);
plot(time(1:end-1), z(2,:), 'LineWidth', 2);
% title('z2 mode (Sliding Left) Over Time');
xlabel('Time (s)');
ylabel('z2');
grid on;

% Third subplot for z3 over time
subplot(3, 1, 3);
plot(time(1:end-1), z(3,:), 'LineWidth', 2);
% title('z3 mode (Sliding Right) Over Time');
xlabel('Time (s)');
ylabel('z3');
grid on;

%% Plots for solver times

figure;
plot(mpc_timestamps, solver_times, 'LineWidth', 2);
hold on;
plot(mpc_timestamps, gurobi_solve_times, 'LineWidth', 2);
xlabel('Simulation Time (s)');
ylabel('Solver Time (s)');
legend('Time By Matlab', 'Time By knitro');
grid on;

%% Path error metrics
% 
% path_error = path_error_MIQP(x, x_star);
% 
% %% plots of path errors
% figure;
% 
% subplot(2,1,1);
% plot(time(1:end), path_error, 'LineWidth', 2);
% % title('Path error x-y');
% xlabel('Time (s)');
% ylabel('$\sqrt{\bar{x}^2 + \bar{y}^2}$ (m)', 'Interpreter', 'latex');
% grid on;
% 
% subplot(2,1,2);
% plot(time(1:end-1), sqrt(dx(1,:).^2 + dx(2,:).^2), 'LineWidth', 2);
% % title('Norm of dx and dy over time');
% xlabel('Time (s)');
% ylabel('$\sqrt{\bar{x}^2 + \bar{y}^2}$ (m)', 'Interpreter', 'latex');
% grid on;

% %% Video area
% % --- Precompute contact points for each time step ---
% contact_x_world = zeros(1, length(x(1,:)));
% contact_y_world = zeros(1, length(x(1,:)));
% 
% for i = 1:length(x(1,:))
%     [x_c, y_c, ~, ~, ~] = calculate_r_c(x(4,i), len, radius, wid, object_shape);
%     R = [cos(x(3,i)), -sin(x(3,i)); sin(x(3,i)), cos(x(3,i))];
%     contact_point_global = R * [x_c; y_c] + [x(1,i); x(2,i)];
% 
%     contact_x_world(i) = contact_point_global(1);
%     contact_y_world(i) = contact_point_global(2);
% end
% 
% % --- 1) Create and configure the video writer
% videoFilename = 'simulation_s_shape.avi';
% video = VideoWriter(videoFilename);
% video.FrameRate = 15;  % Adjust frame rate as desired
% open(video);
% 
% % Create a figure for the animation
% fig = figure();
% 
% % --- 2) Loop over each time step to create frames
% for i = 1:10:length(x(1,:))
% 
%     % Clear figure and hold on for multiple plots
%     clf;  
%     hold on;  
%     grid on;  
%     axis equal;
% 
%     % (Optional) Set the axes to a fixed range if desired
%     % xlim([-1 5]); ylim([-1 5]);  % adjust to your data
% 
%     % --- (A) Plot the completed portion of the object's trajectory so far
%     plot(x(1,1:i), x(2,1:i), 'b-', 'LineWidth', 2, ...
%          'DisplayName','Object Trajectory');
% 
%     % Mark the start and (current) end
%     plot(x(1,1), x(2,1), 'go', 'MarkerFaceColor','g', ...
%          'DisplayName','Start');
%     plot(x(1,i), x(2,i), 'yo', 'MarkerFaceColor','y', ...
%          'DisplayName','Current Position');
% 
%     % --- (B) Plot the desired trajectory up to the current index
%     plot(x_star(1,1:i), x_star(2,1:i), 'k-', 'LineWidth', 2, ...
%          'DisplayName','Desired Trajectory');
% 
%     % Mark the start and end of desired trajectory
%     plot(x_star(1,1),   x_star(2,1),   'go', 'MarkerFaceColor','g');
%     plot(x_star(1,end), x_star(2,end), 'yo', 'MarkerFaceColor','y');
% 
%     % --- (C) Plot the capsule (object) shape at the current step
%     if object_shape == "rectangular_capsule_prism"
%         object = get_capsule_shape(len, radius, x(1,i), x(2,i), x(3,i));
%         fill(object(:,1), object(:,2), 'b', 'FaceAlpha', 0.2, ...
%              'DisplayName','Object Shape');
%     elseif object_shape == "rectangular_prism"
%         object = get_rectangle_shape(len, wid, x(1,i), x(2,i), x(3,i));
%         fill(object(:,1), object(:,2), 'b', 'FaceAlpha', 0.2, ...
%              'DisplayName','Object Shape');
%     end
% 
%     % --- (D) Plot the contact point trajectory up to current index
%     % Assumes you have already computed and stored contact_x_world, contact_y_world
%     % for each time step in arrays of the same length as x.
%     plot(contact_x_world(1:i), contact_y_world(1:i), 'r-', 'LineWidth', 2, ...
%          'DisplayName','Contact Path');
% 
%     % --- (E) Add legend and labels
%     xlabel('x (m)'); ylabel('y (m)');
%     % legend('Location','Best');
% 
%     % --- (F) Capture this frame and write to video
%     drawnow;
%     frame = getframe(fig);
%     writeVideo(video, frame);
% end
% 
% % --- 3) Close the video writer
% close(video);