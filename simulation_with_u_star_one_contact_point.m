%========================================================================%
% Τhis script aims to run a simulation in order to see how the simulation
% runs using u_star as input
%========================================================================%

%% Set up simulation parameters
clear
close all
clc

% len = 0.205;
% wid = 0.155;
% radius = 0.075;
% height = 0.18;
% mass = 4;
len = 0.1;
radius = 0.05;
wid = radius*2;
height = 0.05;
mass = 4;
rectangular_prism_mass = mass * (2*len*radius) / (2*len*radius + pi*radius^2);
cylinder_mass = mass * (pi*radius^2) / (2*len*radius + pi*radius^2);
object_shape = "rectangular_capsule_prism";  % rectangular_capsule_prism or rectangular_prism
I_object = calculate_inertia_matrix(len, wid, height, radius, ...
                                    rectangular_prism_mass, cylinder_mass/2, ...
                                    object_shape);
contact_area = calculate_contact_area(len, wid, radius, object_shape);

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

x_0 = 0.6;
y_0 = -0.2;
x_f = 1;
y_f = 0.2;
v_constant = 0.04;
timestep = 0.002;
trajectory_radius = 0.2;

v_constant_s = 0.055;
x_center = x_0 - trajectory_radius;
y_center = y_0;

% [x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         fifth_trajectory_straight_line(duration, x_0, x_f, y_0, y_f, timestep);

[x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
                        constant_velocity_trajectory_straight_line(v_constant, x_0, x_f, y_0, y_f, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         quarter_circle_trajectory(duration, trajectory_radius, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                     constant_velocity_quarter_circle_trajectory(trajectory_radius, v_constant_s, timestep, x_center, y_center);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                           semi_circle_trajectory(duration, trajectory_radius, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                     constant_velocity_semi_circle_trajectory(trajectory_radius, v_constant_s, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         s_shape_trajectory(duration, trajectory_radius, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                     constant_velocity_s_shape_trajectory(trajectory_radius, v_constant_s, timestep, x_center, y_center);

[fn_star, ft_star, phi_star_dot, phi_star] = ...
                          calculate_u_star_one_contact_point(L, x_star_dot, y_star_dot, theta_star, ...
                          theta_star_dot, len, radius, timestep, duration, wid, object_shape);

x_star = [x_star; y_star; theta_star; phi_star];
u_star = [fn_star; ft_star; phi_star_dot];
x(3,1) =theta_star(1);
x(4,1) = phi_star(1);
x_dot(4,1) = phi_star_dot(1);
%% Run simulation

% Set the control input
dp = [0; 0; 0];
dp_star = [0; 0; 0];
time = 0;

ground_friction_star = zeros(3,1);
ground_friction = zeros(3,1);

for i = 1:duration/timestep

    u = u_star(:,i);

    % Calculate parameters for the motion equation
    w = calculate_motion_model_parameters(u, x(3,i), len, radius, x(4,i), wid, object_shape);
    ground_friction_parameter = 1;
    [gr_frict, number] = calculate_friction_with_ground(L, dp(:,i), ground_friction_parameter);

    if number == 1
        j(i) = 1;
    elseif number == 0
        j(i) = 0;
    end
    
    ground_friction_parameter = 1;
    [gr_fric_star, number] = calculate_friction_with_ground(L, dp_star(:,i), ground_friction_parameter);

    if number == 1
        j_star(i) = 1;
    elseif number == 0
        j_star(i) = 0;
    end

    ground_friction_star(:,i) = -gr_fric_star;

    ground_friction(:,i) = -gr_frict;
    wrench(:,i) = w;

    x_dot(4,i+1) = u(3);
    x(4, i+1) = x(4, i) + x_dot(4, i+1) * timestep;

    x_ddot(1:3, i+1) = inv(diag([mass mass I_object(3,3)])) * (-gr_frict + w);
    x_dot(1:3, i+1) = x_dot(1:3, i) + x_ddot(1:3, i+1) .* timestep;
    x(1:3, i+1) = x(1:3, i) + x_dot(1:3, i+1) .* timestep;

    time(i+1) = time(i) + timestep;
    dp(:,i+1) = [x_dot(1, i+1); x_dot(2, i+1); x_dot(3, i+1)];
    dp_star(:,i+1) = [x_star_dot(i+1); y_star_dot(i+1); theta_star_dot(i+1)];

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
