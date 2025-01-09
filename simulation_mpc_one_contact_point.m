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
mu = 0.5;

L = [1/(mu_ground*F_N)^2 0 0;
     0 1/(mu_ground*F_N)^2 0;
     0 0 1/(alpha*R*mu_ground*F_N)^2];

duration = 8;
x_0 = 0;
x_f = 0.04 * duration;
y_0 = 0;
y_f = 0.04 * duration;

timestep = 0.001;
mpc_timestep = 0.03;
timestep_parameter = mpc_timestep/timestep;
N = 30;
trajectory_radius = 0.2;
v_constant = 0.055;

% [x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         fifth_trajectory_straight_line(duration, x_0, x_f, y_0, y_f, timestep);

% [x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         constant_velocity_trajectory_straight_line(duration, x_0, x_f, y_0, y_f, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         quarter_circle_trajectory(duration, trajectory_radius, timestep);

[x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
                    constant_velocity_quarter_circle_trajectory(trajectory_radius, v_constant, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                           semi_circle_trajectory(duration, trajectory_radius, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                     constant_velocity_semi_circle_trajectory(trajectory_radius, v_constant, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         s_shape_trajectory(duration, trajectory_radius, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                     constant_velocity_s_shape_trajectory(trajectory_radius, v_constant, timestep);

% NOTE: change the initial guess that depends on the trajectory selected
[fn_star, ft_star, phi_star_dot, phi_star, ~] = ...
                          calculate_u_star_one_contact_point(L, x_star_dot, y_star_dot, theta_star, ...
                          theta_star_dot, len, radius, timestep, duration);

x_star = [x_star; y_star; theta_star; phi_star];
u_star = [fn_star; ft_star; phi_star_dot];
% Define the extension for x_star and u_star
x_star_extension = repmat(x_star(:, end), 1, N*timestep_parameter); % Repeat last column of x_star N times
u_star_extension = zeros(size(u_star, 1), N*timestep_parameter);    % Create zero matrix for u_star
% Append the extensions to x_star and u_star
x_star = [x_star, x_star_extension];
u_star = [u_star, u_star_extension];

% System's parameters initialization
x = [0; 0; 0; 0];
x_dot = [0; 0; 0; 0];
x_ddot = [0; 0; 0];
u = [0; 0; 0];
x(:,1) = [trajectory_radius+0.04, 0, 0, phi_star(1)];
% x(:,1) = [0.03, -0.03, 0, phi_star(1)];
x_start = [];

% MPC controller tunable parameters
Q = 60 * diag([4, 4, 0.1, 0]);      % State cost matrix
QN = 26000 * diag([5, 5, 0.1, 0]);   % Terminal state cost matrix
R = 0.01 * diag([1, 1, 0.1]);          % Input cost matrix
%% Run simulation

% Set the control input
dp = [0; 0; 0];
time = 0;

ground_friction = zeros(3,1);
is_start = 1;

for i = 1:floor(duration/timestep)
    dx(:,i) = x(:,i) - x_star(:,i);
    
    if (mod(i*timestep,mpc_timestep) == 0) || i == 1
        [x_star_mpc, u_star_mpc] = create_mpc_star_input(x_star, u_star, ...
                                    N, i, timestep_parameter);
        tic;
        mpc_output = solve_MPC_MIQP(x_star_mpc, u_star_mpc, dx(:,i), mu, L, ...
                                           radius, len, N, mpc_timestep,  Q, QN, R, x_start, is_start);
        solver_times(i) = toc;
        is_start = 0;
        fprintf('Time is: %g\n', time(i));
        fprintf('Error is: %2.4f %2.4f %2.4f %2.4f\n', dx(1,i), dx(2,i), dx(3,i), dx(4,i));
        % fprintf('Velocity is: %2.4f %2.4f %2.4f\n', dp(1,i), dp(2,i), dp(3,i));
        du(:,i) = mpc_output(4*(N+1)+1 : 4*(N+1)+3);
        z(:,i) = mpc_output(4*(N+1)+3*N+1 : 4*(N+1)+3*N+3);
        u(:,i) = du(:,i) + u_star(:,i);
        x_start = mpc_output;
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

    time(i+1) = time(i) + timestep;
    dp(:,i+1) = [x_dot(1, i+1); x_dot(2, i+1); x_dot(3, i+1)];
    

end

%% Present calculation metrics
metrics = calculate_metrics(x, x_star, u, solver_times, N);
disp('Simulation Metrics:');
disp(metrics);

%% Plot Capsule Shape Along Trajectory with Contact Points
figure;
plot(x(1,:), x(2,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Trajectory');
title('Trajectory of x and y Over Time with Capsule Shape and Contact Points');
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

contact_x = [];
contact_y = [];

% Loop to plot contact points and unit vectors separately
for i = 1:length(x(1,:))
    [x_c, y_c, ~, n_c, t_c] = calculate_r_c(x(4,i), len, radius);

    R = [cos(x(3,i)), -sin(x(3,i)); sin(x(3,i)), cos(x(3,i))];
    contact_point_global = R * [x_c; y_c] + [x(1,i); x(2,i)];

    contact_x = [contact_x; contact_point_global(1)];
    contact_y = [contact_y; contact_point_global(2)];

end

% Logical indices for j == 1 and j == 0
idx_1 = (j == 1);
idx_0 = (j == 0);

% Plot for j == 1 (red circles)
plot(contact_x(idx_1), contact_y(idx_1), 'ro', 'MarkerSize', 2);

% Plot for j == 0 (green circles)
plot(contact_x(idx_0), contact_y(idx_0), 'go', 'MarkerSize', 2);

% plot(contact_x, contact_y, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Contact Path');

% Plot x_star and y_star trajectory
plot(x_star(1,:), x_star(2,:), 'k-', 'LineWidth', 2); % x_star trajectory in black
plot(x_star(1,1), x_star(2,1), 'go', 'MarkerFaceColor', 'g'); % Start point of x_star in green
plot(x_star(1,end), x_star(2,end), 'yo', 'MarkerFaceColor', 'y'); % End point of x_star in yellow

legend([capsule_shape_handle; findobj(gca, 'DisplayName', 'Trajectory'); ...
        findobj(gca, 'DisplayName', 'Start'); findobj(gca, 'DisplayName', 'End')], ...
        'Location', 'Best');
        % findobj(gca, 'DisplayName', 'Contact Path')], ...
        % 'Location', 'Best');

hold off;

%% Plot x-x_star, y-y_star and theta-theta_star

figure;

% First subplot for x and x_star over time
subplot(3, 1, 1);
plot(time, x(1, :), 'b-', 'LineWidth', 2); % Plot x
hold on;
plot(time, x_star(1, 1:length(x(1,:))), 'r-', 'LineWidth', 2); % Plot x_star
title('x and x^* Over Time');
xlabel('Time (s)');
ylabel('x (m)');
legend('x', 'x^*', 'Location', 'best');
grid on;

% Second subplot for y and y_star over time
subplot(3, 1, 2);
plot(time, x(2, :), 'b-', 'LineWidth', 2); % Plot y
hold on;
plot(time, x_star(2, 1:length(x(2,:))), 'r-', 'LineWidth', 2); % Plot y_star
title('y and y^* Over Time');
xlabel('Time (s)');
ylabel('y (m)');
legend('y', 'y^*', 'Location', 'best');
grid on;

% Third subplot for theta and theta_star over time
subplot(3, 1, 3);
plot(time, x(3, :), 'b-', 'LineWidth', 2); % Plot theta
hold on;
plot(time, x_star(3, 1:length(x(2,:))), 'r-', 'LineWidth', 2); % Plot theta_star
title('\theta and \theta^* Over Time');
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
title('x\_dot and x\_dot^* Over Time');
xlabel('Time (s)');
ylabel('x\_dot (m/s)');
legend('x\_dot', 'x\_dot^*', 'Location', 'best');
grid on;

% Second subplot for y_dot and y_star_dot over time
subplot(3, 1, 2);
plot(time, x_dot(2,:), 'b-', 'LineWidth', 2); % Plot y_dot
hold on;
plot(time, y_star_dot, 'r-', 'LineWidth', 2); % Plot y_star_dot
title('y\_dot and y\_dot^* Over Time');
xlabel('Time (s)');
ylabel('y\_dot (m/s)');
legend('y\_dot', 'y\_dot^*', 'Location', 'best');
grid on;

% Third subplot for theta_dot and theta_star_dot over time
subplot(3, 1, 3);
plot(time, x_dot(3,:), 'b-', 'LineWidth', 2); % Plot theta_dot
hold on;
plot(time, theta_star_dot, 'r-', 'LineWidth', 2); % Plot theta_star_dot
title('\theta\_dot and \theta\_dot^* Over Time');
xlabel('Time (s)');
ylabel('\theta\_dot (rad/s)');
legend('\theta\_dot', '\theta\_dot^*', 'Location', 'best');
grid on;

%% Plot dx dy and dtheta in seperate plots

figure;
% First subplot for dx over time
subplot(3, 1, 1);
plot(time(1:end-1), dx(1,:), 'LineWidth', 2);
title('dx Over Time');
xlabel('Time (s)');
ylabel('dx (m)');
grid on;

% Second subplot for dy over time
subplot(3, 1, 2);
plot(time(1:end-1), dx(2,:), 'LineWidth', 2);
title('dy Over Time');
xlabel('Time (s)');
ylabel('dy (m)');
grid on;

% Third subplot for dtheta over time
subplot(3, 1, 3);
plot(time(1:end-1), dx(3,:), 'LineWidth', 2);
title('dtheta Over Time');
xlabel('Time (s)');
ylabel('dtheta (rad)');
grid on;

%% Plot fn ft and phi_dot in seperate plots

figure;
% First subplot for fn over time
subplot(3, 1, 1);
plot(time(1:end-1), u(1,:), 'LineWidth', 2);
title('fn Over Time');
xlabel('Time (s)');
ylabel('fn (N)');
grid on;

% Second subplot for ft over time
subplot(3, 1, 2);
plot(time(1:end-1), u(2,:), 'LineWidth', 2);
title('ft Over Time');
xlabel('Time (s)');
ylabel('ft (N)');
grid on;

% Third subplot for phi_dot over time
subplot(3, 1, 3);
plot(time(1:end-1), u(3,:), 'LineWidth', 2);
title('phi dot Over Time');
xlabel('Time (s)');
ylabel('phi dot (rad/sec)');
grid on;

%% Plot dfn dft and dphi_dot in seperate plots

figure;
% First subplot for dfn over time
subplot(3, 1, 1);
plot(time(1:length(du(1,:))), du(1,:), 'LineWidth', 2);
title('dfn Over Time');
xlabel('Time (s)');
ylabel('dfn (N)');
grid on;

% Second subplot for dft over time
subplot(3, 1, 2);
plot(time(1:length(du(2,:))), du(2,:), 'LineWidth', 2);
title('dft Over Time');
xlabel('Time (s)');
ylabel('dft (N)');
grid on;

% Third subplot for dphi_dot over time
subplot(3, 1, 3);
plot(time(1:length(du(3,:))), du(3,:), 'LineWidth', 2);
title('dphi dot Over Time');
xlabel('Time (s)');
ylabel('dphi dot (rad/sec)');
grid on;

%% Plot z1 z2 and z3 in seperate plots

figure;
% First subplot for z1 over time
subplot(3, 1, 1);
plot(time(1:end-1), z(1,:), 'LineWidth', 2);
title('z1 mode (sticking) Over Time');
xlabel('Time (s)');
ylabel('binary z1');
grid on;

% Second subplot for z2 over time
subplot(3, 1, 2);
plot(time(1:end-1), z(2,:), 'LineWidth', 2);
title('z2 mode (Sliding Left) Over Time');
xlabel('Time (s)');
ylabel('binary z2');
grid on;

% Third subplot for z3 over time
subplot(3, 1, 3);
plot(time(1:end-1), z(3,:), 'LineWidth', 2);
title('z3 mode (Sliding Right) Over Time');
xlabel('Time (s)');
ylabel('binary z3');
grid on;

%% Plot phi and phi_dot in seperate plots

figure;
% First subplot for dfn over time
subplot(2, 1, 1);
plot(time(1:end), x(4,:), 'LineWidth', 2);
title('phi Over Time');
xlabel('Time (s)');
ylabel('phi (rad)');
grid on;

% Second subplot for phi_dot over time
subplot(2, 1, 2);
plot(time(1:end), x_dot(4,:), 'LineWidth', 2);
title('phi dot Over Time');
xlabel('Time (s)');
ylabel('phi dot (rad/sec)');
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

% Fourth subplot for torque on z-axis
subplot(4,1,4);
plot(time(1:end-1), sqrt(wrench(1,:).^2 + wrench(2,:).^2 + wrench(3,:).^2), 'LineWidth', 2);
hold on;
plot(time(1:end-1), sqrt(ground_friction(1,:).^2 + ground_friction(2,:).^2 + ground_friction(3,:).^2), 'LineWidth', 2);
title('Norm of friction force over time');
xlabel('Time (s)');
ylabel('Torque (Nm)');
legend('Wrench', 'Ground Friction');
grid on;

sgtitle('Wrench and Ground Friction Forces and Torque over Time');

%% Plots for debug

figure;
plot(time(1:length(time)-1), j, 'LineWidth', 2);
title('binary debug');
xlabel('Time (s)');
ylabel('on-off mode');
legend('Ground Friction');
grid on;
