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
mass = 4;
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
mu = 0.3;

L = [1/(mu_ground*F_N)^2 0 0;
     0 1/(mu_ground*F_N)^2 0;
     0 0 1/(alpha*R*mu_ground*F_N)^2];

duration = 5;
x_0 = 0;
x_f = 0.04 * duration;
y_0 = 0;
y_f = 0.04 * duration;

timestep = 0.001;
mpc_timestep = 0.04;
timestep_parameter = mpc_timestep/timestep;
control_frequency = 0.04;
N = 20; 
trajectory_radius = 0.2;
v_constant = 0.055;
x_center = 0;
y_center = 0;

% [x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         fifth_trajectory_straight_line(duration, x_0, x_f, y_0, y_f, timestep);

% [x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         constant_velocity_trajectory_straight_line(duration, x_0, x_f, y_0, y_f, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         quarter_circle_trajectory(duration, trajectory_radius, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                     constant_velocity_quarter_circle_trajectory(trajectory_radius, v_constant, timestep, x_center, y_center);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                           semi_circle_trajectory(duration, trajectory_radius, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                     constant_velocity_semi_circle_trajectory(trajectory_radius, v_constant, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         s_shape_trajectory(duration, trajectory_radius, timestep);

[x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
                    constant_velocity_s_shape_trajectory(trajectory_radius, v_constant, timestep, x_center, y_center);

% NOTE: change the initial guess that depends on the trajectory selected
[fn_star, ft_star, phi_star_dot, phi_star, ~] = ...
                          calculate_u_star_one_contact_point(L, x_star_dot, y_star_dot, theta_star, ...
                          theta_star_dot, len, radius, timestep, duration);

x_star = [x_star; y_star; theta_star; phi_star];
u_star = [fn_star; ft_star; phi_star_dot];
% u_star = [fn_star - 2 * ones(size(fn_star)); ft_star; phi_star_dot];
% Define the extension for x_star and u_star
x_star_extension = repmat(x_star(:, end), 1, control_frequency/timestep + N*timestep_parameter); % Repeat last column of x_star N times
u_star_extension = zeros(size(u_star, 1), control_frequency/timestep + N*timestep_parameter);    % Create zero matrix for u_star
% u_star_extension = repmat(u_star(:, end), 1, control_frequency/timestep + N*timestep_parameter);
% Append the extensions to x_star and u_star
x_star = [x_star, x_star_extension];
u_star = [u_star, u_star_extension];
% u_star(1, :) = 22 * ones(size(x_star(1, :)));
% u_star(2:3, :) = zeros(size(x_star(2:3, :)));
% u_star = zeros(size(x_star(1:3, :)));

% System's parameters initialization
x = [0; 0; 0; 0];
x_dot = [0; 0; 0; 0];
x_ddot = [0; 0; 0];
u = [0; 0; 0];
x(:,1) = [trajectory_radius+0.04, -0.04, 0, 3*pi/2];
% x(:,1) = [0.06, -0.06, 0, 3*pi/2];
[x_c, y_c, ~, n_c, t_c] = calculate_r_c(x(4,1), len, radius);
x_pc_world(:,1) = x(1:2, 1) + [cos(x(3, 1)) -sin(x(3, 1)); sin(x(3, 1)) cos(x(3, 1))] *[x_c; y_c];
mpc_output = [];

% MPC controller tunable parameters for 0.04s
% Q = 80 * diag([5, 5, 0.1, 0]);
% QN = 42000 * diag([5, 5, 0.1, 0]);
% R = 0.02 * diag([1, 1, 0.01]);

% MPC controller tunable parameters for 0.04s
Q = 100 * diag([5, 5, 0.01, 0]);
QN = 40000 * diag([5, 5, 0.01, 0]);
R = 0.04 * diag([1, 1, 0.01]);

%% Run simulation

% Set the control input
dp = [0; 0; 0];
time = 0;
mpc_timestamps = 0;
solver_times = 0;

ground_friction = zeros(3,1);
is_start = 1;
k = 1;

for i = 1:floor(duration/timestep)
    dx(:,i) = x(:,i) - x_star(:,i);
    
    if i == 1
        simulation_type_flag = true;
        [x_star_mpc, u_star_mpc, dx_mpc] = create_mpc_star_constant_u(x_star, u_star,...
                                    N, i, timestep_parameter, control_frequency,...
                                    u(:,i), x(:,i), len, radius, dp(:,i), timestep, ...
                                    L, mass, I_object, simulation_type_flag);
        tic;
        [mpc_output, gurobi_solve_time] = solve_MPC_MIQP(x_star_mpc, u_star_mpc, dx(:,i), mu, L, ...
                                           radius, len, N, mpc_timestep,  Q, QN, R, mpc_output, is_start);
        solver_times(k) = round(toc*1000) / 1000;
        gurobi_solve_times(k) = round(gurobi_solve_time*1000) / 1000;
        mpc_timestamps(k) = time(i);
        % is_start = 0;
        fprintf('Time is: %g\n', time(i));
        fprintf('Error is: %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f\n', dx(1,i), dx(2,i), dx(3,i), dx(4,i), ... 
                     sqrt(dx(1,i)^2 + dx(2,i)^2), sqrt(dx(1,i)^2 + dx(2,i)^2 + dx(3,i)^2));
        du(:,i) = mpc_output(4*(N+1)+1 : 4*(N+1)+3);
        z(:,i) = mpc_output(4*(N+1)+3*N+1 : 4*(N+1)+3*N+3);
        u(:,i) = du(:,i) + u_star(:,i);
        k = k + 1;
    elseif (mod(i*timestep, control_frequency) == 0 || i == 2)
        [x_star_mpc, u_star_mpc, dx_mpc] = create_mpc_star_constant_u(x_star, u_star,...
                                    N, i, timestep_parameter, control_frequency,...
                                    u(:,i-1), x(:,i), len, radius, dp(:,i), timestep, ...
                                    L, mass, I_object, simulation_type_flag);
        du(:,i) = mpc_output(4*(N+1)+1 : 4*(N+1)+3);
        z(:,i) = mpc_output(4*(N+1)+3*N+1 : 4*(N+1)+3*N+3);
        u(:,i) = du(:,i) + u_star(:,i);
        tic;
        [mpc_output, gurobi_solve_time] = solve_MPC_MIQP(x_star_mpc, u_star_mpc, dx_mpc, mu, L, ...
                                           radius, len, N, mpc_timestep,  Q, QN, R, mpc_output, is_start);
        solver_times(k) = round(toc*1000) / 1000;
        gurobi_solve_times(k) = round(gurobi_solve_time*1000) / 1000;
        mpc_timestamps(k) = time(i);
        fprintf('Time is: %g\n', time(i));
        fprintf('Error is: %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f\n', dx(1,i), dx(2,i), dx(3,i), dx(4,i), ... 
                     sqrt(dx(1,i)^2 + dx(2,i)^2), sqrt(dx(1,i)^2 + dx(2,i)^2 + dx(3,i)^2));
        k = k + 1;
    else
        u(:,i) = u(:,i-1);
        z(:,i) = z(:,i-1);
    end
    
    % Calculate parameters for the motion equation
    w = calculate_motion_model_parameters(u(:,i), x(3,i), len, radius, x(4,i));
    ground_friction_parameter = 1;
    [gr_frict, number] = calculate_friction_with_ground(L, dp(:,i), ground_friction_parameter);

    if number == 1
        j(i) = 1;
    elseif number == 0
        j(i) = 0;
    end

    ground_friction(:,i) = -gr_frict;
    wrench(:,i) = w;

    x_dot(4,i+1) = u(3,i);
    x(4, i+1) = x(4, i) + x_dot(4, i+1) * timestep;

    % x_ddot(1:3, i+1) = inv(diag([mass mass I_object(3,3)])) * (-gr_frict + w);
    x_ddot(1:3, i+1) = diag([mass mass I_object(3,3)]) \ (-gr_frict + w);
    x_dot(1:3, i+1) = x_dot(1:3, i) + x_ddot(1:3, i+1) .* timestep;
    x(1:3, i+1) = x(1:3, i) + x_dot(1:3, i+1) .* timestep;

    [v_pc_body(:,i), v_pc_world(:,i)] = calculate_robot_velocity(L, x(:, i), len, radius, u(:,i));
    x_pc_world(:,i+1) = x_pc_world(:,i) + v_pc_world(:,i)*timestep;

    time(i+1) = time(i) + timestep;
    dp(:,i+1) = [x_dot(1, i+1); x_dot(2, i+1); x_dot(3, i+1)];
    
end

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

capsule_shape_handle = [];

% Plot the object shape at several points along the trajectory
for i = 1:round(length(x(1,:))/10):length(x(1,:))
    capsule_shape = get_capsule_shape(len, radius, x(1,i), x(2,i), x(3,i));

    if isempty(capsule_shape_handle) 
        capsule_shape_handle = fill(capsule_shape(:,1), capsule_shape(:,2), 'b', 'FaceAlpha', 0.2, 'DisplayName', 'Capsule Shape');
    else
        fill(capsule_shape(:,1), capsule_shape(:,2), 'b', 'FaceAlpha', 0.2);
    end
end

contact_x_world = [];
contact_y_world = [];

% Loop to plot contact points and unit vectors separately
for i = 1:length(x(1,:))
    [x_c, y_c, ~, n_c, t_c] = calculate_r_c(x(4,i), len, radius);

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
plot(contact_x_world(idx_0), contact_y_world(idx_0), 'ro', 'MarkerSize', 2);

% plot(contact_x, contact_y, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Contact Path');

% Plot x_star and y_star trajectory
plot(x_star(1,:), x_star(2,:), 'k-', 'LineWidth', 2, 'DisplayName', 'Desired trajectory'); % x_star trajectory in black
plot(x_star(1,1), x_star(2,1), 'go', 'MarkerFaceColor', 'g'); % Start point of x_star in green
plot(x_star(1,end), x_star(2,end), 'yo', 'MarkerFaceColor', 'y'); % End point of x_star in yellow

legend([capsule_shape_handle; findobj(gca, 'DisplayName', 'Trajectory'); ...
        findobj(gca, 'DisplayName', 'Desired trajectory'); ...
        findobj(gca, 'DisplayName', 'Start'); findobj(gca, 'DisplayName', 'End'); ...
        contact_handle], 'Location', 'Best');

hold off;

%% Plot x-x_star, y-y_star and theta-theta_star

figure;

% First subplot for x and x_star over time
subplot(3, 1, 1);
plot(time, x(1, :), 'b-', 'LineWidth', 2); % Plot x
hold on;
plot(time, x_star(1, 1:length(x(1,:))), 'r-', 'LineWidth', 2); % Plot x_star
% title('x and x^* Over Time');
xlabel('Time (s)');
ylabel('x (m)');
legend('x', 'x^*', 'Location', 'best');
grid on;

% Second subplot for y and y_star over time
subplot(3, 1, 2);
plot(time, x(2, :), 'b-', 'LineWidth', 2); % Plot y
hold on;
plot(time, x_star(2, 1:length(x(2,:))), 'r-', 'LineWidth', 2); % Plot y_star
% title('y and y^* Over Time');
xlabel('Time (s)');
ylabel('y (m)');
legend('y', 'y^*', 'Location', 'best');
grid on;

% Third subplot for theta and theta_star over time
subplot(3, 1, 3);
plot(time, x(3, :), 'b-', 'LineWidth', 2); % Plot theta
hold on;
plot(time, x_star(3, 1:length(x(2,:))), 'r-', 'LineWidth', 2); % Plot theta_star
% title('\theta and \theta^* Over Time');
xlabel('Time (s)');
ylabel('\theta (rad)');
legend('\theta', '\theta^*', 'Location', 'best');
grid on;

%% Plot x_dot-x_star_dot, y_dot-y_star_dot, and theta_dot-theta_star_dot

figure;

% First subplot for x_dot and x_star_dot over time
subplot(3, 1, 1);
plot(time, x_dot(1,:), 'b-', 'LineWidth', 2); % Plot x_dot
hold on;
plot(time, x_star_dot(1,:), 'r-', 'LineWidth', 2); % Plot x_star_dot
xlabel('Time (s)');
ylabel('$\dot{x}$ (m/s)', 'Interpreter', 'latex');
legend({'$\dot{x}$', '$\dot{x}^*$'}, 'Interpreter', 'latex', 'Location', 'best');
grid on;

% Second subplot for y_dot and y_star_dot over time
subplot(3, 1, 2);
plot(time, x_dot(2,:), 'b-', 'LineWidth', 2); % Plot y_dot
hold on;
plot(time, y_star_dot, 'r-', 'LineWidth', 2); % Plot y_star_dot
xlabel('Time (s)');
ylabel('$\dot{y}$ (m/s)', 'Interpreter', 'latex');
legend({'$\dot{y}$', '$\dot{y}^*$'}, 'Interpreter', 'latex', 'Location', 'best');
grid on;

% Third subplot for theta_dot and theta_star_dot over time
subplot(3, 1, 3);
plot(time, x_dot(3,:), 'b-', 'LineWidth', 2); % Plot theta_dot
hold on;
plot(time, theta_star_dot, 'r-', 'LineWidth', 2); % Plot theta_star_dot
xlabel('Time (s)');
ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter', 'latex');
legend({'$\dot{\theta}$', '$\dot{\theta}^*$'}, 'Interpreter', 'latex', 'Location', 'best');
grid on;


%% Plot dx dy and dtheta in seperate plots

figure;
% First subplot for dx over time
subplot(3, 1, 1);
plot(time(1:end-1), dx(1,:), 'LineWidth', 2);
% title('dx Over Time');
xlabel('Time (s)');
ylabel('$\bar{x}$ (m)', 'Interpreter', 'latex');
grid on;

% Second subplot for dy over time
subplot(3, 1, 2);
plot(time(1:end-1), dx(2,:), 'LineWidth', 2);
% title('dy Over Time');
xlabel('Time (s)');
ylabel('$\bar{y}$ (m)', 'Interpreter', 'latex');
grid on;

% Third subplot for dtheta over time
subplot(3, 1, 3);
plot(time(1:end-1), dx(3,:), 'LineWidth', 2);
% title('dtheta Over Time');
xlabel('Time (s)');
ylabel('$\bar{\theta}$ (rad)', 'Interpreter', 'latex');
grid on;

figure;
subplot(2,1,1);
plot(time(1:end-1), sqrt(dx(1,:).^2 + dx(2,:).^2), 'LineWidth', 2);
% title('Norm of dx and dy over time');
xlabel('Time (s)');
ylabel('$\sqrt{\bar{x}^2 + \bar{y}^2}$ (m)', 'Interpreter', 'latex');
grid on;

subplot(2,1,2);
plot(time(1:end-1), sqrt(dx(1,:).^2 + dx(2,:).^2 + dx(3,:).^2), 'LineWidth', 2);
% title('Norm of dx and dy dtheta over time');
xlabel('Time (s)');
ylabel('$\sqrt{\bar{x}^2 + \bar{y}^2 + \bar{\theta}^2}$', 'Interpreter', 'latex');
grid on;

%% Plot fn ft and phi_dot in seperate plots

figure;
% First subplot for fn over time
subplot(3, 1, 1);
plot(time(1:end-1), u(1,:), 'LineWidth', 2);
% title('fn Over Time');
xlabel('Time (s)');
ylabel('fn (N)');
grid on;

% Second subplot for ft over time
subplot(3, 1, 2);
plot(time(1:end-1), u(2,:), 'LineWidth', 2);
% title('ft Over Time');
xlabel('Time (s)');
ylabel('ft (N)');
grid on;

% Third subplot for phi_dot over time
subplot(3, 1, 3);
plot(time(1:end-1), u(3,:), 'LineWidth', 2);
% title('phi dot Over Time');
xlabel('Time (s)');
ylabel('$\dot{\phi}$ (rad/sec)', 'Interpreter', 'latex');
grid on;


%% Plot dfn dft and dphi_dot in seperate plots

figure;

% First subplot for dfn over time
subplot(3, 1, 1);
plot(time(1:length(du(1,:))), du(1,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\bar{f}_n$ (N)', 'Interpreter', 'latex');
grid on;

% Second subplot for dft over time
subplot(3, 1, 2);
plot(time(1:length(du(2,:))), du(2,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\bar{f}_t$ (N)', 'Interpreter', 'latex');
grid on;

% Third subplot for dphi_dot over time
subplot(3, 1, 3);
plot(time(1:length(du(3,:))), du(3,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{\bar{\phi}}$ (rad/s)', 'Interpreter', 'latex');
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

%% Plot phi and phi_dot in seperate plots

figure;

% First subplot for phi over time
subplot(2, 1, 1);
plot(time(1:end), x(4,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\bar{\phi}$ (rad)', 'Interpreter', 'latex');
grid on;

% Second subplot for phi_dot over time
subplot(2, 1, 2);
plot(time(1:end), x_dot(4,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{\bar{\phi}}$ (rad/s)', 'Interpreter', 'latex');
grid on;


%% Plot wrench and ground friction over time over time

figure;
% First subplot for force on x-axis
subplot(4,1,1);
plot(time(1:end-1), wrench(1,:), 'LineWidth', 2);
hold on;
plot(time(1:end-1), ground_friction(1,:), 'LineWidth', 2);
title('Force on x-axis over time');
xlabel('Time (s)');
ylabel('Force (N)');
legend('Wrench', 'Ground Friction');
grid on;

% Second subplot for force on y-axis
subplot(4,1,2);
plot(time(1:end-1), wrench(2,:), 'LineWidth', 2);
hold on;
plot(time(1:end-1), ground_friction(2,:), 'LineWidth', 2);
title('Force on y-axis over time');
xlabel('Time (s)');
ylabel('Force (N)');
legend('Wrench', 'Ground Friction');
grid on;

% Third subplot for torque on z-axis
subplot(4,1,3);
plot(time(1:end-1), wrench(3,:), 'LineWidth', 2);
hold on;
plot(time(1:end-1), ground_friction(3,:), 'LineWidth', 2);
title('Torque on z-axis over time');
xlabel('Time (s)');
ylabel('Torque (Nm)');
legend('Wrench', 'Ground Friction');
grid on;

% Fourth subplot for norm of force
subplot(4,1,4);
plot(time(1:end-1), sqrt(wrench(1,:).^2 + wrench(2,:).^2 + wrench(3,:).^2), 'LineWidth', 2);
hold on;
plot(time(1:end-1), sqrt(ground_friction(1,:).^2 + ground_friction(2,:).^2 + ground_friction(3,:).^2), 'LineWidth', 2);
title('Norm of friction force over time');
xlabel('Time (s)');
ylabel('Torque (Nm)');
legend('Wrench', 'Ground Friction');
grid on;

%% Plots for debug

figure;
plot(time(1:length(time)-1), j, 'LineWidth', 2);
title('binary debug');
xlabel('Time (s)');
ylabel('on-off mode');
grid on;

%% Plots for solver times

figure;
plot(mpc_timestamps, solver_times, 'LineWidth', 2);
hold on;
plot(mpc_timestamps, gurobi_solve_times, 'LineWidth', 2);
% title('Solver time');
xlabel('Simulation Time (s)');
ylabel('Solver Time (s)');
legend('Time By Matlab', 'Time By Gurobi');
grid on;

%% Plot the velocity of the robitic arm world frame

figure;
% First subplot for velocity on x-axis world frame and object frame
subplot(3,1,1);
plot(time(1:end-1), v_pc_world(1,:), 'LineWidth', 2);
hold on;
plot(time(1:end-1), contact_x_dot_world(2:end), 'LineWidth', 2);
title('Velocity of robot on x-axis world frame');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('Velocity x-axis robot', 'Velocity x-axis contact point');
grid on;

% Second subplot for velocity on y-axis world frame and object frame
subplot(3,1,2);
plot(time(1:end-1), v_pc_world(2,:), 'LineWidth', 2);
hold on;
plot(time(1:end-1), contact_y_dot_world(2:end), 'LineWidth', 2);
title('Velocity of robot on y-axis world frame');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('Velocity y-axis robot', 'Velocity y-axis contact point');
grid on;

% Second subplot for norm velocity world frame and object frame
subplot(3,1,3);
plot(time(1:end-1), sqrt(v_pc_world(1,:).^2 + v_pc_world(2,:).^2), 'LineWidth', 2);
hold on;
plot(time(1:end-1), sqrt(contact_x_dot_world(2:end).^2 + contact_y_dot_world(2:end).^2), 'LineWidth', 2);
title('Norm of velocity of robot and contact point world frame');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('Norm y-axis robot', 'Norm y-axis contact point');
grid on;

%% Plot the velocity of the robitic arm body frame

figure;
% First subplot for velocity on x-axis world frame and object frame
subplot(3,1,1);
plot(time(1:end-1), v_pc_body(1,:), 'LineWidth', 2);
hold on;
plot(time(1:end-1), contact_x_dot_body(2:end), 'LineWidth', 2);
title('Velocity of robot on x-axis body frame');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('Velocity x-axis robot', 'Velocity x-axis contact point');
grid on;

% Second subplot for velocity on y-axis world frame and object frame
subplot(3,1,2);
plot(time(1:end-1), v_pc_body(2,:), 'LineWidth', 2);
hold on;
plot(time(1:end-1), contact_y_dot_body(2:end), 'LineWidth', 2);
title('Velocity of robot on y-axis body frame');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('Velocity y-axis robot', 'Velocity y-axis contact point');
grid on;

% Second subplot for norm velocity world frame and object frame
subplot(3,1,3);
plot(time(1:end-1), sqrt(v_pc_body(1,:).^2 + v_pc_body(2,:).^2), 'LineWidth', 2);
hold on;
plot(time(1:end-1), sqrt(contact_x_dot_body(2:end).^2 + contact_y_dot_body(2:end).^2), 'LineWidth', 2);
title('Norm of velocity of robot and contact point body frame');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('Norm y-axis robot', 'Norm y-axis contact point');
grid on;

%% Plot the position of the robitic arm over time

figure;
% First subplot for potistion on x-axis world frame
subplot(3,1,1);
plot(time(1:end), x_pc_world(1,:), 'LineWidth', 2);
hold on;
plot(time(1:end), contact_x_world, 'LineWidth', 2);
title('Position of robot on x-axis over time');
xlabel('Time (s)');
ylabel('Position (m)');
legend('positions contact point x-axis', 'positions object x-axis');
grid on;

% Second subplot for potistion on y-axis world frame
subplot(3,1,2);
plot(time(1:end), x_pc_world(2,:), 'LineWidth', 2);
hold on;
plot(time(1:end), contact_y_world, 'LineWidth', 2);
title('Position of robot on y-axis over time');
xlabel('Time (s)');
ylabel('Position (m)');
legend('positions contact point y-axis', 'positions object y-axis');
grid on;

% Third subplot for norm of potistion world frame
subplot(3,1,3);
plot(x_pc_world(1,:), x_pc_world(2,:), 'LineWidth', 2);
hold on;
plot(contact_x_world, contact_y_world, 'LineWidth', 2);
title('2D plot of position of robot over time');
xlabel('x_c (m)');
ylabel('y_c (m)');
grid on;

%% Path error metrics

% [path_error, x_y_error, theta_error] = path_error_MIQP(x, x_star);
% 
% %% plots of path errors
% figure;
% % First subplot for path error
% subplot(3,1,1);
% plot(time(1:end), path_error, 'LineWidth', 2);
% % title('Path error x-y-theta');
% xlabel('Time (s)');
% ylabel('$\sqrt{\bar{x}^2 + \bar{y}^2 + \bar{\theta}^2}$', 'Interpreter', 'latex');
% grid on;
% 
% % Second subplot for path error
% subplot(3,1,2);
% plot(time(1:end), x_y_error, 'LineWidth', 2);
% % title('Path error x-y');
% xlabel('Time (s)');
% ylabel('$\sqrt{\bar{x}^2 + \bar{y}^2}$ (m)', 'Interpreter', 'latex');
% grid on;
% 
% % Third subplot for path error
% subplot(3,1,3);
% plot(time(1:end), theta_error, 'LineWidth', 2);
% % title('Path error theta');
% xlabel('Time (s)');
% ylabel('$\bar{\theta}$ (rad)', 'Interpreter', 'latex');
% 
% grid on;
% 
% 

