%========================================================================%
% This script calculates u_star based on a specific trajectory and creates 
% many plots that gives a good understanding of the solutions
%========================================================================%

%% Set up simulation parameters
clear
close all
clc

len = 0.205;
wid = 0.155;
radius = 0.075;
height = 0.18;
mass = 4;
% len = 0.1;
% radius = 0.05;
% wid = radius*2;
% height = 0.05;
% mass = 4;
rectangular_prism_mass = mass * (2*len*radius) / (2*len*radius + pi*radius^2);
cylinder_mass = mass * (pi*radius^2) / (2*len*radius + pi*radius^2);
object_shape = "rectangular_prism";  % rectangular_capsule_prism or rectangular_prism
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
y_0 = -0.3;
x_f = 1;
y_f = 0.2;
v_constant = 0.04;
timestep = 0.002;
trajectory_radius = 0.2;

v_constant_s = 0.055;
x_center = x_0 - trajectory_radius;
y_center = y_0;

mpc_timestep = 0.04;
timestep_parameter = mpc_timestep/timestep;
control_frequency = 0.04;
N = 20; 

% [x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, time] = ...
%                         fifth_trajectory_straight_line(duration, x_0, x_f, y_0, y_f, timestep);

[x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, time, duration] = ...
                        constant_velocity_trajectory_straight_line(v_constant, x_0, x_f, y_0, y_f, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, time] = ...
%                         quarter_circle_trajectory(duration, trajectory_radius, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, time, duration] = ...
%                     constant_velocity_quarter_circle_trajectory(trajectory_radius, v_constant_s, timestep, x_center, y_center);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, time] = ...
%                           semi_circle_trajectory(duration, trajectory_radius, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, time, duration] = ...
%                     constant_velocity_semi_circle_trajectory(trajectory_radius, v_constant_s, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, time] = ...
%                         s_shape_trajectory(duration, trajectory_radius, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, time, duration] = ...
%                     constant_velocity_s_shape_trajectory(trajectory_radius, v_constant_s, timestep, x_center, y_center);

[fn_star, ft_star, phi_star_dot, phi_star] = ...
                          calculate_u_star_one_contact_point(L, x_star_dot, y_star_dot, theta_star, ...
                          theta_star_dot, len, radius, timestep, duration, wid, object_shape);

for i = 1:length(phi_star_dot)
    if phi_star_dot(i) ~= 0
        phi_star_dot(i) = 0;
    end
end

T = table(x_star(:), y_star(:), theta_star(:), phi_star(:), fn_star(:), ft_star(:), phi_star_dot(:), ...
          'VariableNames', {'x_star', 'y_star', 'theta_star', 'phi_star', 'fn_star', 'ft_star', 'phi_star_dot'});

writetable(T, 'simulation_data_straight_line.csv');
% writetable(T, 'simulation_data_quarter_circle_20cm_radius.csv');
% writetable(T, 'simulation_data_semi_circle_20cm_radius.csv');
% writetable(T, 'simulation_data_s_shape_20cm_radius.csv')


%% plot x_c and y_c

for i = 1:(duration/timestep)+1
    [x_c, y_c, r_c, n_c, t_c] = calculate_r_c(phi_star(i), len, radius, wid, object_shape);
    plot_y_c(i) = y_c;
    plot_x_c(i) = x_c;
end

figure;
% First subplot for y_c over time
subplot(2, 1, 1);
plot(time, plot_y_c, 'Linewidth', 2);
title('y_c Over Time');
xlabel('Time (s)');
ylabel('y_c (m)');
grid on;

% Second subplot for x_c over time
subplot(2, 1, 2);
plot(time, plot_x_c, 'Linewidth', 2);
title('x_c Over Time');
xlabel('Time (s)');
ylabel('x_c (m)');
grid on;

%% plot x_star y_star theta_star and the whole 2d surface

figure;
plot(x_star, y_star, 'Linewidth', 2);
title('Trajectory of x and y Over Time');
xlabel('x (m)');
ylabel('y (m)');
grid on;
axis equal;

hold on;
plot(x_star(1), y_star(1), 'go', 'MarkerFaceColor', 'g');
plot(x_star(end), y_star(end), 'yo', 'MarkerFaceColor', 'y');
legend('Trajectory', 'Start', 'End');
hold off;

figure;
% First subplot for x_star over time
subplot(3, 1, 1);
plot(time, x_star, 'Linewidth', 2);
title('x star Over Time');
xlabel('Time (s)');
ylabel('x star (m)');
grid on;

% Second subplot for y_star over time
subplot(3, 1, 2);
plot(time, y_star, 'Linewidth', 2);
title('y star Over Time');
xlabel('Time (s)');
ylabel('y star (m)');
grid on;

% Third subplot for theta_star over time
subplot(3, 1, 3);
plot(time, theta_star, 'Linewidth', 2);
title('theta star Over Time');
xlabel('Time (s)');
ylabel('theta star (rad)');
grid on;

%% Plot x_star and y_star and theta_star velocity over time

figure;
% First subplot for x_star_dot over time
subplot(3, 1, 1);
plot(time, x_star_dot, 'Linewidth', 2);
title('x star Velocity Over Time');
xlabel('Time (s)');
ylabel('x star dot (m/s)');
grid on;

% Second subplot for y_star_dot over time
subplot(3, 1, 2);
plot(time, y_star_dot, 'Linewidth', 2);
title('y star Velocity Over Time');
xlabel('Time (s)');
ylabel('y star dot (m/s)');
grid on;

% Third subplot for theta_star_dot over time
subplot(3, 1, 3);
plot(time, theta_star_dot, 'Linewidth', 2);
title('theta star Velocity Over Time');
xlabel('Time (s)');
ylabel('theta star dot(rad/s)');
grid on;

%% Plot phi_ phi_star

figure;
% First subplot for phi_star over time
subplot(2, 1, 1);
plot(time, phi_star, 'Linewidth', 2);
title('phi star Over Time');
xlabel('Time (s)');
ylabel('phi star (rad)');
grid on;

% Second subplot for phi_star_dot over time
subplot(2, 1, 2);
plot(time, phi_star_dot, 'Linewidth', 2);
title('phi star Velocity Over Time');
xlabel('Time (s)');
ylabel('phi star dot (rad/s)');
grid on;


%% Plot fn_star and ft_star

figure;
% First subplot for fn_star over time
subplot(2, 1, 1);
plot(time, fn_star, 'Linewidth', 2);
title('fn star Over Time');
xlabel('t (s)');
ylabel('fn (N)');
grid on;

% Second subplot for ft_star over time
subplot(2, 1, 2);
plot(time, ft_star, 'Linewidth', 2);
title('ft star Over Time');
xlabel('t (s)');
ylabel('ft (N)');
grid on;

%% Plot Capsule Shape Along Trajectory with Contact Points
figure;
plot(x_star, y_star, 'Linewidth', 2, 'DisplayName', 'Trajectory');
title('Trajectory of x and y Over Time with Capsule Shape and Contact Points');
xlabel('x (m)');
ylabel('y (m)');
grid on;
axis equal;

hold on;

plot(x_star(1), y_star(1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(x_star(end), y_star(end), 'yo', 'MarkerFaceColor', 'y', 'DisplayName', 'End');

capsule_shape_handle = [];

% Plot the object shape at several points along the trajectory
for i = 1:round(length(x_star)/10):length(x_star)
    capsule_shape = get_capsule_shape(len, radius, x_star(i), y_star(i), theta_star(i));

    if isempty(capsule_shape_handle) 
        capsule_shape_handle = fill(capsule_shape(:,1), capsule_shape(:,2), 'b', 'FaceAlpha', 0.3, 'DisplayName', 'Capsule Shape');
    else
        fill(capsule_shape(:,1), capsule_shape(:,2), 'b', 'FaceAlpha', 0.3);
    end
end

contact_x = [];
contact_y = [];

% Define a scaling factor for the vectors to make them visible
vector_scale = 0.01;

% Loop to plot contact points and unit vectors separately
for i = 1:length(x_star)
    [x_c, y_c, ~, n_c, t_c] = calculate_r_c(phi_star(i), len, radius, wid, object_shape);

    R = [cos(theta_star(i)), -sin(theta_star(i)); sin(theta_star(i)), cos(theta_star(i))];
    contact_point_global = R * [x_c; y_c] + [x_star(i); y_star(i)];

    contact_x = [contact_x; contact_point_global(1)];
    contact_y = [contact_y; contact_point_global(2)];

    % % Transform the unit vectors to global coordinates
    % n_c_global = R * n_c;
    % t_c_global = R * t_c;
    % 
    % % Plot the normal and tangent unit vectors at the contact point
    % quiver(contact_point_global(1), contact_point_global(2), ...
    %        n_c_global(1) * vector_scale, n_c_global(2) * vector_scale, ...
    %        'g', 'Linewidth', 1.5, 'MaxHeadSize', 1, 'DisplayName', 'Normal Vector');
    % 
    % quiver(contact_point_global(1), contact_point_global(2), ...
    %        t_c_global(1) * vector_scale, t_c_global(2) * vector_scale, ...
    %        'm', 'Linewidth', 1.5, 'MaxHeadSize', 1, 'DisplayName', 'Tangent Vector');
end

plot(contact_x, contact_y, 'r-', 'Linewidth', 1.5, 'DisplayName', 'Contact Path');

legend([capsule_shape_handle; findobj(gca, 'DisplayName', 'Trajectory'); ...
        findobj(gca, 'DisplayName', 'Start'); findobj(gca, 'DisplayName', 'End'); ...
        findobj(gca, 'DisplayName', 'Contact Path')], ...
        'Location', 'Best');

hold off;

%% Plot the torque created by fn star and ft star

for i = 1:length(fn_star)
    [x_c, y_c, ~, n_c, t_c] = calculate_r_c(phi_star(i), len, radius, wid, object_shape);
    torque_fn(i) = -[y_c -x_c] * n_c .* fn_star(i);
    torque_ft(i) = -[y_c -x_c] * t_c .* ft_star(i);
end

figure;
plot(time, torque_fn, 'Linewidth', 2);
title('fn torque star Over Time');
xlabel('t (s)');
ylabel('fn torque (Nm)');
grid on;

figure;
plot(time, torque_ft, 'Linewidth', 2);
title('ft torque star Over Time');
xlabel('t (s)');
ylabel('ft torque (Nm)');
grid on;






