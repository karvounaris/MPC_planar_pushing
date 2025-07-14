%========================================================================%
% Τhis script runs a simulation using the mpc controller
%========================================================================%

%% Set up simulation parameters
clear
close all
clc

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

v_constant_s = 0.055;
x_center = x_0 - trajectory_radius;
y_center = y_0;

mpc_timestep = 0.04;
N = 30; 

% [x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                         constant_velocity_trajectory_straight_line(v_constant, x_0, x_f, y_0, y_f, timestep);

[x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
                    constant_velocity_quarter_circle_trajectory(trajectory_radius, v_constant_s, timestep, x_center, y_center);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                     constant_velocity_semi_circle_trajectory(trajectory_radius, v_constant_s, timestep, x_center, y_center);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                     constant_velocity_s_shape_trajectory(trajectory_radius, v_constant_s, timestep, x_center, y_center);

x_star = [x_star; y_star; theta_star; ones(size(theta_star))*3*pi/2];
% Define the extension for x_star and u_star
x_star_extension = repmat(x_star(:, end), 1, round(2*N*mpc_timestep/timestep));
% Append the extensions to x_star and u_star
x_star = [x_star, x_star_extension];

% System's parameters initialization
x = [0; 0; 0; 0];
x_dot = [0; 0; 0; 0];
x_ddot = [0; 0; 0];
u = [0; 0; 0; 0];
% x(:,1) = [x_0 + 0.05, y_0 - 0.05, 0, 5*pi/4];
% x(:,1) = [x_0 + 0.05, y_0 - 0.05, 0, 7*pi/4];
% x(:,1) = [x_0 + 0.05, y_0 - 0.05, 0, 5*pi/4];
% x(:,1) = [x_0 + 0.05, y_0 - 0.05, 0, 7*pi/4];
x(:,1) = [x_0 + 0.05, y_0 - 0.05, 0, 3*pi/2];
% obstacle_data = [[0.1; 0.2; 0.05], [0.5; 0.1; 0.05], [0.3; 0.4; 0.05]];  % [x_obstacle, y_obstacle, radius_obstacle]
% obstacle_data = [-0.12; 0.04; 0.04];  % [x_obstacle, y_obstacle, radius_obstacle]
obstacle_data = [];
n_obstacles = size(obstacle_data, 2);

% for the straight line change the upper and bottom limits for phi_dot to -3 +3
% MPC controller tunable parameters
Q = 60 * diag([4, 4, 0.0001, 0]);
QN = 8000 * diag([4, 4, 0.005, 0]);
R = 0.0001 * diag([1, 1, 0, 0]);

%% Run simulation
% profile on
% Set the control input
dp = [0; 0; 0];
time = zeros(1, floor(duration/timestep)+1);
time(1) = 0;
mpc_timestamps = 0;
solver_times = 0;
knitro_solve_times = 0;
last_time_MIQP_start = 0;
w_eps0 = 50;
k_eps = 0.2;

X_init = repmat(x(:,1), 1, N+1);
U_init = ones(4, N);
E_init = zeros(1, N);

[solver_controller, args] = create_MPC_MPCC_controller( ...
    Q, QN, R, w_eps0, k_eps, N, mpc_timestep, mu, L, radius, len, wid, object_shape, n_obstacles);

args.x0 = [reshape(X_init, 4*(N+1),1); reshape(U_init, 4*N,1); E_init'];

ground_friction = zeros(3,1);
is_start = 1;
k = 1;

for i = 1:floor(duration/timestep)

    dx(:,i) = x(:,i) - x_star(:,i);
    
    if i == 1
        last_time_MIQP_start = i * timestep;
        for j = 1:N+1
            x_star_mpc(:,j) = x_star(:, i + round((j-1)*mpc_timestep/timestep));
        end
        args.p = [x(:,i); x_star_mpc(:); obstacle_data(:)];
        tic;
        sol = solver_controller('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx, ...
                     'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p);
        solver_times(k) = round(toc*1000) / 1000;
        info = solver_controller.stats();
        knitro_solve_times(k) = round(info.t_wall_total * 1000) / 1000;
        mpc_timestamps(k) = time(i);
        fprintf('Time is: %g\n', time(i));
        fprintf('Error is: %2.4f %2.4f %2.4f %2.4f %2.4f\n', dx(1,i), dx(2,i), dx(3,i), dx(4,i), ... 
                     sqrt(dx(1,i)^2 + dx(2,i)^2));
        X_opt = reshape(sol.x(1:4*(N+1)), 4, N+1);
        U_opt = reshape(sol.x(4*(N+1)+1:4*(N+1)+4*N), 4, N);
        Eps_opt = reshape(sol.x(4*(N+1)+4*N+1:end), 1, N);
        X_init = full([X_opt(:,2:end), X_opt(:,end)]);
        U_init = full([U_opt(:,2:end), U_opt(:,end)]);
        E_init = full([Eps_opt(2:end), Eps_opt(end)]);
        args.x0 = [reshape(X_init, 4*(N+1),1); reshape(U_init, 4*N,1); E_init'];
        eps(i) = full(Eps_opt(1));
        u(:,i) = full(U_opt(:,1));
        k = k + 1;
    elseif (i * timestep - last_time_MIQP_start > solver_times(end) || i == 2)
        last_time_MIQP_start = i * timestep;
        for j = 1:N+1
            x_star_mpc(:,j) = x_star(:, i + round((j-1)*mpc_timestep/timestep));
        end
        args.p = [x(:,i); x_star_mpc(:); obstacle_data(:)]; 
        eps(i) = full(Eps_opt(1));
        u(:,i) = full(U_opt(:,1));
        tic;
        sol = solver_controller('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx, ...
                     'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p);
        solver_times(k) = round(toc*1000) / 1000;
        info = solver_controller.stats();
        knitro_solve_times(k) = round(info.t_wall_total * 1000) / 1000;
        mpc_timestamps(k) = time(i);
        fprintf('Time is: %g\n', time(i));
        fprintf('Error is: %2.4f %2.4f %2.4f %2.4f %2.4f\n', dx(1,i), dx(2,i), dx(3,i), dx(4,i), ... 
                     sqrt(dx(1,i)^2 + dx(2,i)^2));
        if info.unified_return_status == "SOLVER_RET_SUCCESS"
            X_opt = reshape(sol.x(1:4*(N+1)), 4, N+1);
            U_opt = reshape(sol.x(4*(N+1)+1:4*(N+1)+4*N), 4, N);
            Eps_opt = reshape(sol.x(4*(N+1)+4*N+1: end), 1, N);
            X_init = full([X_opt(:,2:end), X_opt(:,end)]);
            U_init = full([U_opt(:,2:end), U_opt(:,end)]);
            E_init = full([Eps_opt(2:end), Eps_opt(end)]);
            args.x0 = [reshape(X_init, 4*(N+1),1); reshape(U_init, 4*N,1); E_init'];
        else
            X_init = [X_init(:,2:end), X_init(:,end)];
            U_init = [U_init(:,2:end), U_init(:,end)];
            E_init = [E_init(2:end), E_init(end)];
            args.x0 = [reshape(X_init, 4*(N+1),1); reshape(U_init, 4*N,1); E_init'];
        end
        k = k + 1;
    elseif i * timestep - last_time_MIQP_start > mpc_timestep
        X_opt = X_init;
        U_opt = U_init;
        Eps_opt = E_init;
        X_init = full([X_opt(:,2:end), X_opt(:,end)]);
        U_init = full([U_opt(:,2:end), U_opt(:,end)]);
        E_init = full([Eps_opt(2:end), Eps_opt(end)]);
        args.x0 = [reshape(X_init, 4*(N+1),1); reshape(U_init, 4*N,1); E_init'];
        eps(i) = full(Eps_opt(1));
        u(:,i) = full(U_opt(:,1));
    else
        u(:,i) = u(:,i-1);
        eps(:,i) = eps(i-1);
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

    x_dot(4,i+1) = u(3,i)-u(4,i);
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

contact_handle = plot(contact_x_world, contact_y_world, 'r-', 'LineWidth', 2, 'DisplayName', 'Contact Point');

for i = 1:round(length(x(1,:))/5):length(x(1,:))

    if object_shape == "rectangular_capsule_prism"
        shape_handle = get_capsule_shape(len, radius, x(1,i), x(2,i), x(3,i));
    elseif object_shape == "rectangular_prism"
        shape_handle = get_rectangle_shape(len, radius, x(1,i), x(2,i), x(3,i));
    end

    if i == 1 
        shape_handle_handle = fill(shape_handle(:,1), shape_handle(:,2), 'b', 'FaceAlpha', 0.2, 'DisplayName', 'Capsule Shape');
    else
        fill(shape_handle(:,1), shape_handle(:,2), 'b', 'FaceAlpha', 0.2);
    end

    r_contact = 0.01;
    theta_circ = linspace(0, 2*pi, 50);
    [x_c, y_c, ~, n_c, t_c] = calculate_r_c(x(4,i), len, radius, wid, object_shape);
    R = [cos(x(3,i)), -sin(x(3,i)); sin(x(3,i)), cos(x(3,i))];
    contact_point_global_graph = R * [x_c; y_c-r_contact] + [x(1,i); x(2,i)];

    contact_x_world_graph = [];
    contact_y_world_graph = [];

    contact_x_world_graph = [contact_x_world_graph; contact_point_global_graph(1)];
    contact_y_world_graph = [contact_y_world_graph; contact_point_global_graph(2)];
    
    circle_x = r_contact * cos(theta_circ) + contact_x_world_graph(end);
    circle_y = r_contact * sin(theta_circ) + contact_y_world_graph(end);
    
    if i == 1 
        shape_handle_handle_cycle = fill(circle_x, circle_y, 'r', 'FaceAlpha', 1.0, 'EdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', 'Robot');
    else
        fill(circle_x, circle_y, 'r', 'FaceAlpha', 1.0, 'EdgeColor', 'k', 'LineWidth', 0.5);
    end
end

% Plot x_star and y_star trajectory
plot(x_star(1,:), x_star(2,:), 'k-', 'LineWidth', 2, 'DisplayName', 'Desired trajectory'); % x_star trajectory in black
plot(x_star(1,1), x_star(2,1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(x_star(1,end), x_star(2,end), 'yo', 'MarkerFaceColor', 'y', 'DisplayName', 'End');

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

proxy_contact_point = plot(nan, nan, 'ro', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', ...
                           'LineWidth', 0.5, 'DisplayName', 'Robot');

if n_obstacles >=1
    legend([shape_handle_handle; ...
        proxy_contact_point; ...
        findobj(gca, 'DisplayName', 'Trajectory'); ...
        findobj(gca, 'DisplayName', 'Desired trajectory'); ...
        findobj(gca, 'MarkerFaceColor', 'g'); ...
        findobj(gca, 'MarkerFaceColor', 'y'); ...
        contact_handle, obstacle_handles(1)], ...
        'Location', 'Best');
else
    legend([shape_handle_handle; ...
        proxy_contact_point; ...
        findobj(gca, 'DisplayName', 'Trajectory'); ...
        findobj(gca, 'DisplayName', 'Desired trajectory'); ...
        findobj(gca, 'MarkerFaceColor', 'g'); ...
        findobj(gca, 'MarkerFaceColor', 'y'); ...
        contact_handle], ...
        'Location', 'Best');
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
subplot(4, 1, 1);
plot(time(1:end-1), u(1,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('fn (N)');
grid on;

subplot(4, 1, 2);
plot(time(1:end-1), u(2,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('ft (N)');
grid on;

subplot(4, 1, 3);
plot(time(1:end-1), u(3,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{\phi_+}$ (rad/sec)', 'Interpreter', 'latex');
grid on;

subplot(4, 1, 4);
plot(time(1:end-1), u(4,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{\phi_-}$ (rad/sec)', 'Interpreter', 'latex');
grid on;

%% Plots for epsilon

figure;
plot(time(1:end-1), eps, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('epsilon value');
grid on;

%% Plots for solver times

figure;
plot(mpc_timestamps, solver_times, 'LineWidth', 2);
hold on;
plot(mpc_timestamps, knitro_solve_times, 'LineWidth', 2);
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