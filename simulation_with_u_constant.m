%========================================================================%
% Τhis script aims to run a simulation in order to test the environment. We
% will not use control algorithms, instead we will use a constant input u 
% to the system and make it work as expected from the physics perspective.
%========================================================================%

%% Set up simulation parameters
clear
close all
clc

% Object's parameters
len = 0.06;
radius = 0.03;
width = radius*2;
height = 0.03;
mass = 1;
rectangular_prism_mass = mass * (2*len*radius) / (2*len*radius + pi*radius^2);
cylinder_mass = mass * (pi*radius^2) / (2*len*radius + pi*radius^2);
object_shape = "rectangular_capsule_prism";
I_object = calculate_inertia_matrix(len, width, height, radius, ...
                                    rectangular_prism_mass, cylinder_mass/2, ...
                                    object_shape);
contact_area = calculate_contact_area(len, width, radius, object_shape);

% Contact points' parameters
C = 1;
% phi_C_start = [11*pi/8, 13*pi/8];
% phi_C_start = 11*pi/8;
phi_C_start = 13*pi/8;
% phi_C_start = 3*pi/2;
% phi_C_start = 0;
contact_points = cell(C,4);
for i = 1:C
    contact_points{i,1} = 0;            % x_C
    contact_points{i,2} = 0;            % y_C
    contact_points{i,3} = 0;            % phi_C
    contact_points{i,4} = zeros(2,3);   % J_C
    contact_points{i,5} = [0; 0];       % n_C
    contact_points{i,6} = [0; 0];       % t_C
end
N = zeros(3,C);
T = zeros(3,C);

% System's parameters initialization
x = 1;
x_dot = 0;
x_double_dot = 0;
y = 1;
y_dot = 0;
y_double_dot = 0;
theta = 0;
theta_dot = 0;
theta_double_dot = 0;
phi_C = zeros(C,1);
for i = 1:C
    phi_C(i,1) = phi_C_start(i);
end
phi_C_dot = zeros(C,1);
mu_ground = 0.5;          % coefficient of friction ground_object


%% Run simulation

% Simulation parameters
duration = 10;
timestep = 0.03;

% Set the control input
u = [zeros(C,1); zeros(C,1); zeros(C,1)];
dp = [0; 0; 0];
time = 0;

for i = 1:duration/timestep + 1

    % Calculate parameters for the motion equation
    [R, contact_points, w, ground_friction] = calculate_motion_model_parameters(u, theta(i), len, radius, C, ...
                                                contact_points, N, T, mass, mu_ground, contact_area, dp, phi_C_start);

    for j = 1:C
        phi_C_dot(j, i+1) = u(2*C+j,1);
        phi_C(j, i+1) = phi_C(j, i) + phi_C_dot(j, i+1) * timestep;
    end

    x_double_dot(i+1) = (-ground_friction(1) + w(1))/mass;
    % x_double_dot(i+1) = w(1)/mass;
    x_dot(i+1) = x_dot(i) + x_double_dot(i+1) * timestep;
    x(i+1) = x(i) + x_dot(i+1) * timestep;

    y_double_dot(i+1) = (-ground_friction(2) + w(2))/mass;
    % y_double_dot(i+1) = w(2)/mass;
    y_dot(i+1) = y_dot(i) + y_double_dot(i+1) * timestep;
    y(i+1) = y(i) + y_dot(i+1) * timestep;

    theta_double_dot(i+1) = (-ground_friction(3) + w(3))/I_object(3,3);
    % theta_double_dot(i+1) = w(3)/I_object(3,3);
    theta_dot(i+1) = theta_dot(i) + theta_double_dot(i+1) * timestep;
    theta(i+1) = theta(i) + theta_dot(i+1) * timestep;
    
    u = [0.3 * ones(C,1); zeros(C,1); zeros(C,1)];    
    time(i+1) = time(i) + timestep;
    dp = [x(i+1) - x(i); y(i+1) - y(i); theta(i+1) - theta(i)];

end

%% Plot results

% Plot x and y position together in a 2D grid
figure;
plot(x, y, 'LineWidth', 2);
title('Trajectory of x and y Over Time');
xlabel('x (m)');
ylabel('y (m)');
grid on;

hold on;
plot(x(1), y(1), 'go', 'MarkerFaceColor', 'g'); % Start point in green
plot(x(end), y(end), 'yo', 'MarkerFaceColor', 'y'); % End point in red
legend('Trajectory', 'Start', 'End');
hold off;

% % Plot x velocity over time
% figure;
% plot(time, x_dot, 'LineWidth', 2);
% title('x Velocity Over Time');
% xlabel('Time (s)');
% ylabel('x Velocity (m/s)');
% grid on;
% 
% % Plot y velocity over time
% figure;
% plot(time, y_dot, 'LineWidth', 2);
% title('y Velocity Over Time');
% xlabel('Time (s)');
% ylabel('y Velocity (m/s)');
% grid on;
% 
% % Plot x acceleration over time
% figure;
% plot(time, x_double_dot, 'LineWidth', 2);
% title('x Acceleration Over Time');
% xlabel('Time (s)');
% ylabel('x Acceleration (m/s^2)');
% grid on;
% 
% % Plot y acceleration over time
% figure;
% plot(time, y_double_dot, 'LineWidth', 2);
% title('y Acceleration Over Time');
% xlabel('Time (s)');
% ylabel('y Acceleration (m/s^2)');
% grid on;
% Define object shape for plotting (rectangular capsule prism)

%% Plot Capsule Shape Along Trajectory with Contact Points
figure;
plot(x_star, y_star, 'LineWidth', 2, 'DisplayName', 'Trajectory');
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
for i = 1:round(length(x_star)/20):length(x_star)
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
    [x_c, y_c, ~, n_c, t_c] = calculate_r_c(phi_star(i), len, radius);

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
    %        'g', 'LineWidth', 1.5, 'MaxHeadSize', 1, 'DisplayName', 'Normal Vector');
    % 
    % quiver(contact_point_global(1), contact_point_global(2), ...
    %        t_c_global(1) * vector_scale, t_c_global(2) * vector_scale, ...
    %        'm', 'LineWidth', 1.5, 'MaxHeadSize', 1, 'DisplayName', 'Tangent Vector');
end

plot(contact_x, contact_y, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Contact Path');

legend([capsule_shape_handle; findobj(gca, 'DisplayName', 'Trajectory'); ...
        findobj(gca, 'DisplayName', 'Start'); findobj(gca, 'DisplayName', 'End'); ...
        findobj(gca, 'DisplayName', 'Contact Path')], ...
        'Location', 'Best');

hold off;

