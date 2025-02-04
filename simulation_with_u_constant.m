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

% System's parameters initialization
x = [1; 1; pi/4; 3*pi/2];
x_dot = [0; 0; 0; 0];
x_ddot = [0; 0; 0; 0];
mu_ground = 0.5;          % coefficient of friction ground_object


%% Run simulation

% Simulation parameters
duration = 5;
timestep = 0.001;

% Set the control input
u = [20; 0; 0];
dp = [0; 0; 0];
time = 0;

for i = 1:duration/timestep + 1

    % Calculate parameters for the motion equation
    w = calculate_motion_model_parameters(u, x(3,i), len, radius, x(4,i));
    ground_friction_parameter = 1;
    [gr_frict, number] = calculate_friction_with_ground(L, dp(:,i), ground_friction_parameter);

    ground_friction(:,i) = -gr_frict;
    wrench(:,i) = w;

    x_dot(4,i+1) = u(3);
    x(4, i+1) = x(4, i) + x_dot(4, i+1) * timestep;

    % x_ddot(1:3, i+1) = inv(diag([mass mass I_object(3,3)])) * (-gr_frict + w);
    x_ddot(1:3, i+1) = diag([mass mass I_object(3,3)]) \ (-gr_frict + w);
    x_dot(1:3, i+1) = x_dot(1:3, i) + x_ddot(1:3, i+1) .* timestep;
    x(1:3, i+1) = x(1:3, i) + x_dot(1:3, i+1) .* timestep;
    
    time(i+1) = time(i) + timestep;
    dp(:,i+1) = [x_dot(1, i+1); x_dot(2, i+1); x_dot(3, i+1)];    % Calculate parameters for the motion equation

end

%% Plot results

% Plot x and y position together in a 2D grid
figure;
plot(x(1,:), x(2,:), 'LineWidth', 2);
title('Trajectory of x and y Over Time');
xlabel('x (m)');
ylabel('y (m)');
grid on;

hold on;
plot(x(1,1), x(2,1), 'go', 'MarkerFaceColor', 'g'); % Start point in green
plot(x(1,end), x(2,end), 'yo', 'MarkerFaceColor', 'y'); % End point in red
legend('Trajectory', 'Start', 'End');
hold off;

%% Plot x, y and theta

figure;

% First subplot for x over time
subplot(3, 1, 1);
plot(time, x(1, :), 'b-', 'LineWidth', 2); % Plot x
title('x Over Time');
xlabel('Time (s)');
ylabel('x (m)');
grid on;

% Second subplot for y over time
subplot(3, 1, 2);
plot(time, x(2, :), 'b-', 'LineWidth', 2); % Plot y
title('y and Over Time');
xlabel('Time (s)');
ylabel('y (m)');
grid on;

% Third subplot for theta over time
subplot(3, 1, 3);
plot(time, x(3, :), 'b-', 'LineWidth', 2); % Plot theta
title('\theta and Over Time');
xlabel('Time (s)');
ylabel('\theta (rad)');
grid on;

%% Plot x_dot, y_dot, and theta_dot

figure;

% First subplot for x_dot over time
subplot(3, 1, 1);
plot(time, x_dot(1,:), 'b-', 'LineWidth', 2); % Plot x_dot
title('x\_dot and Over Time');
xlabel('Time (s)');
ylabel('x\_dot (m/s)');
grid on;

% Second subplot for y_dot over time
subplot(3, 1, 2);
plot(time, x_dot(2,:), 'b-', 'LineWidth', 2); % Plot y_dot
title('y\_dot and Over Time');
xlabel('Time (s)');
ylabel('y\_dot (m/s)');
grid on;

% Third subplot for theta_dot over time
subplot(3, 1, 3);
plot(time, x_dot(3,:), 'b-', 'LineWidth', 2); % Plot theta_dot
title('\theta\_dot and Over Time');
xlabel('Time (s)');
ylabel('\theta\_dot (rad/s)');
grid on;

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
plot(contact_x_world(idx_1), contact_y_world(idx_1), 'ro', 'MarkerSize', 2);

% Plot for j == 0 (green circles)
plot(contact_x_world(idx_0), contact_y_world(idx_0), 'go', 'MarkerSize', 2);

legend([capsule_shape_handle; findobj(gca, 'DisplayName', 'Trajectory'); ...
        findobj(gca, 'DisplayName', 'Start'); findobj(gca, 'DisplayName', 'End')], ...
        'Location', 'Best');

hold off;