%========================================================================%
% Τhis script aims to run a simulation in order to see how the simulation
% runs using u_star as input
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
x = [0; 0; 0; 0];
x_dot = [0; 0; 0; 0];
x_ddot = [0; 0; 0];

trajectory_radius = 0.3;

% Simulation parameters
duration = 6;
timestep = 0.001;

x_0 = 0;
x_f = 0.04 * duration;
y_0 = 0;
y_f = 0.04 * duration;

% [x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         fifth_trajectory_straight_line(duration, x_0, x_f, y_0, y_f, timestep);

[x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, phi_star, phi_star_dot, ~] = ...
                        constant_velocity_trajectory_straight_line(duration, x_0, x_f, y_0, y_f, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         quarter_circle_trajectory(duration, trajectory_radius, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, phi_star, phi_star_dot, ~, duration] = ...
%                     constant_velocity_quarter_circle_trajectory(trajectory_radius, v_constant, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                           semi_circle_trajectory(duration, trajectory_radius, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                     constant_velocity_semi_circle_trajectory(trajectory_radius, v_constant, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~] = ...
%                         s_shape_trajectory(duration, trajectory_radius, timestep);

% [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, ~, duration] = ...
%                     constant_velocity_s_shape_trajectory(trajectory_radius, v_constant, timestep);

[fn_star, ft_star, phi_star_dot, phi_star, ~] = ...
                          calculate_u_star_one_contact_point(L, x_star_dot, y_star_dot, theta_star, ...
                          theta_star_dot, len, radius, timestep, duration);

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
    w = calculate_motion_model_parameters(u, x(3,i), len, radius, x(4,i));
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

%% Plot results

% Plot x and y position together in a 2D grid
figure;
plot(x(1,:), x(2,:), 'LineWidth', 2);
title('Trajectory of x and y Over Time');
xlabel('x (m)');
ylabel('y (m)');
axis equal;
grid on;

hold on;
plot(x(1,1), x(2,1), 'go', 'MarkerFaceColor', 'g'); % Start point in green
plot(x(1,end), x(2,end), 'yo', 'MarkerFaceColor', 'y'); % End point in red
legend('Trajectory', 'Start', 'End');
hold off;


figure;
% First subplot for x over time
subplot(3, 1, 1);
plot(time, x(1,:), 'LineWidth', 2);
title('x Over Time');
xlabel('Time (s)');
ylabel('x (m)');
grid on;

% Second subplot for y over time
subplot(3, 1, 2);
plot(time, x(2,:), 'LineWidth', 2);
title('y Over Time');
xlabel('Time (s)');
ylabel('y (m)');
grid on;

% Third subplot for theta over time
subplot(3, 1, 3);
plot(time, x(3,:), 'LineWidth', 2);
title('theta Over Time');
xlabel('Time (s)');
ylabel('theta (rad)');
grid on;

%% Plot x and y and theta velocity over time

figure;
% First subplot for x_dot over time
subplot(3, 1, 1);
plot(time, x_dot(1,:), 'LineWidth', 2);
title('x Velocity Over Time');
xlabel('Time (s)');
ylabel('x dot (m/s)');
grid on;

% Second subplot for y_dot over time
subplot(3, 1, 2);
plot(time, x_dot(2,:), 'LineWidth', 2);
title('y Velocity Over Time');
xlabel('Time (s)');
ylabel('y dot (m/s)');
grid on;

% Third subplot for theta_dot over time
subplot(3, 1, 3);
plot(time, x_dot(3,:), 'LineWidth', 2);
title('theta Velocity Over Time');
xlabel('Time (s)');
ylabel('theta dot (rad/s)');
grid on;

%% Plot x_double_dot and y_double_dot and theta_double_dot acceleration over time

figure;
% First subplot for x_double_dot over time
subplot(3, 1, 1);
plot(time, x_ddot(1,:), 'LineWidth', 2);
title('x Acceleration Over Time');
xlabel('Time (s)');
ylabel('x double dot (m/s^2)');
grid on;

% Second subplot for y_double_dot over time
subplot(3, 1, 2);
plot(time, x_ddot(2,:), 'LineWidth', 2);
title('y Acceleration Over Time');
xlabel('Time (s)');
ylabel('y double dot (m/s^2)');
grid on;

% Third subplot for theta_double_dot over time
subplot(3, 1, 3);
plot(time, x_ddot(3,:), 'LineWidth', 2);
title('theta Acceleration Over Time');
xlabel('Time (s)');
ylabel('theta double dot (rad/s^2)');
grid on;

%% Plot phi phi_dot

figure;
% First subplot for phi over time
subplot(2, 1, 1);
plot(time, x(4,:), 'LineWidth', 2);
title('phi Over Time');
xlabel('Time (s)');
ylabel('phi (rad)');
grid on;

% Second subplot for phi_dot over time
subplot(2, 1, 2);
plot(time, x_dot(4,:), 'LineWidth', 2);
title('phi Velocity Over Time');
xlabel('Time (s)');
ylabel('phi dot (rad/s)');
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
        capsule_shape_handle = fill(capsule_shape(:,1), capsule_shape(:,2), 'b', 'FaceAlpha', 0.3, 'DisplayName', 'Capsule Shape');
    else
        fill(capsule_shape(:,1), capsule_shape(:,2), 'b', 'FaceAlpha', 0.3);
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

plot(contact_x, contact_y, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Contact Path');


% Plot x_star and y_star trajectory
plot(x_star(1,:), y_star(1,:), 'k-', 'LineWidth', 2); % x_star trajectory in black
plot(x_star(1,1), y_star(1,1), 'go', 'MarkerFaceColor', 'g'); % Start point of x_star in green
plot(x_star(1,end), y_star(1,end), 'yo', 'MarkerFaceColor', 'y'); % End point of x_star in yellow

legend([capsule_shape_handle; findobj(gca, 'DisplayName', 'Trajectory'); ...
        findobj(gca, 'DisplayName', 'Start'); findobj(gca, 'DisplayName', 'End'); ...
        findobj(gca, 'DisplayName', 'Contact Path')], ...
        'Location', 'Best');

hold off;


%% Plot wrench and ground friction over time over time

figure;
% First subplot for force on x-axis
subplot(3,1,1);
plot(time(1:end-1), wrench(1,:), 'LineWidth', 2);
hold on;
plot(time(1:end-1), ground_friction(1,:), 'LineWidth', 2);
title('Force on x-axis over time');
xlabel('Time (s)');
ylabel('Force (N)');
legend('Wrench', 'Ground Friction');
grid on;

% Second subplot for force on y-axis
subplot(3,1,2);
plot(time(1:end-1), wrench(2,:), 'LineWidth', 2);
hold on;
plot(time(1:end-1), ground_friction(2,:), 'LineWidth', 2);
title('Force on y-axis over time');
xlabel('Time (s)');
ylabel('Force (N)');
legend('Wrench', 'Ground Friction');
grid on;

% Third subplot for torque on z-axis
subplot(3,1,3);
plot(time(1:end-1), wrench(3,:), 'LineWidth', 2);
hold on;
plot(time(1:end-1), ground_friction(3,:), 'LineWidth', 2);
title('Torque on z-axis over time');
xlabel('Time (s)');
ylabel('Torque (Nm)');
legend('Wrench', 'Ground Friction');
grid on;

sgtitle('Wrench and Ground Friction Forces and Torque over Time');


%% Plot ground friction over time over time for simulated and star

figure;

% First subplot for force on x-axis
subplot(3,1,1);
plot(time(1:end-1), ground_friction(1,:), 'LineWidth', 2);
hold on;
plot(time(1:end-1), ground_friction_star(1,:), 'LineWidth', 2);
title('Force on x-axis over time');
xlabel('Time (s)');
ylabel('Force (N)');
legend('Ground Friction', 'Ground Friction Star');
grid on;

% Second subplot for force on y-axis
subplot(3,1,2);
plot(time(1:end-1), ground_friction(2,:), 'LineWidth', 2);
hold on;
plot(time(1:end-1), ground_friction_star(2,:), 'LineWidth', 2);
title('Force on y-axis over time');
xlabel('Time (s)');
ylabel('Force (N)');
legend('Ground Friction', 'Ground Friction Star');
grid on;

% Third subplot for torque on z-axis
subplot(3,1,3);
plot(time(1:end-1), ground_friction(3,:), 'LineWidth', 2);
hold on;
plot(time(1:end-1), ground_friction_star(3,:), 'LineWidth', 2);
title('Torque on z-axis over time');
xlabel('Time (s)');
ylabel('Torque (Nm)');
legend('Ground Friction', 'Ground Friction Star');
grid on;

sgtitle('Ground Friction from star trajectory and simulated trajectory');


%% Plots for debug

figure;
plot(time(1:length(time)-1), j, 'LineWidth', 2);
hold on;
plot(time(1:length(time)-1), j_star, 'LineWidth', 2);
title('binary debug');
xlabel('Time (s)');
ylabel('on-off mode');
legend('Ground Friction', 'Ground Friction Star');
grid on;

%% velocities beween simulated and star
figure;

% First subplot for x velocities
subplot(3,1,1);
plot(time, dp(1,:), 'LineWidth', 2);
hold on;
plot(time, dp_star(1,:), 'LineWidth', 2);
title('Velocity x and Velocity Star x');
xlabel('Time (s)');
ylabel('Velocity x (m/s)');
legend('Velocity x', 'Velocity Star x');
grid on;

% Second subplot for y velocities
subplot(3,1,2);
plot(time, dp(2,:), 'LineWidth', 2);
hold on;
plot(time, dp_star(2,:), 'LineWidth', 2);
title('Velocity y and Velocity Star y');
xlabel('Time (s)');
ylabel('Velocity y (m/s)');
legend('Velocity y', 'Velocity Star y');
grid on;

% Third subplot for theta velocities
subplot(3,1,3);
plot(time, dp(3,:), 'LineWidth', 2);
hold on;
plot(time, dp_star(3,:), 'LineWidth', 2);
title('Velocity theta and Velocity Star theta');
xlabel('Time (s)');
ylabel('Angular Velocity theta (rad/s)');
legend('Angular Velocity theta', 'Angular Velocity Star theta');
grid on;

sgtitle('Velocities between Simulated and Star');

%% paronomastis kai arithmitis
figure;
paronomastis = zeros(size(time));
staronomastis = zeros(size(time));


for i = 1:length(paronomastis)
    paronomastis(i) = sqrt(dp(:, i)'*inv(L)*dp(:, i));
    staronomastis(i) = sqrt(dp_star(:, i)'*inv(L)*dp_star(:, i));
    % paronomastis(i) = dp(:, i)'*inv(L)*dp(:, i);
    % staronomastis(i) = dp_star(:, i)'*inv(L)*dp_star(:, i);
end

plot(time, paronomastis, 'LineWidth', 2);
hold on;
plot(time, staronomastis, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('paronomastis');
legend('Paronomastis', 'Paronomastis Star');
grid on;

arithmitis = zeros(3, length(time));
starithmitis = zeros(3, length(time));


for i = 1:length(arithmitis)
    arithmitis(:, i) = inv(L)*dp(:, i);
    starithmitis(:, i) = inv(L)*dp_star(:, i);
end

figure;

subplot(3,1,1);
plot(time, arithmitis(1,:), 'LineWidth', 2);
hold on;
plot(time, starithmitis(1,:), 'LineWidth', 2);
title('Arithmitis');
xlabel('Time (s)');
ylabel('arithmitis x (N)');
legend('arithmitis x', 'starithmitis x');
grid on;

subplot(3,1,2);
plot(time, arithmitis(2,:), 'LineWidth', 2);
hold on;
plot(time, starithmitis(2,:), 'LineWidth', 2);
title('Arithmitis');
xlabel('Time (s)');
ylabel('trivi y (N)');
legend('arithmitis y', 'starithmitis y');
grid on;

subplot(3,1,3);
plot(time, arithmitis(3,:), 'LineWidth', 2);
hold on;
plot(time, starithmitis(3,:), 'LineWidth', 2);
title('Arithmitis');
xlabel('Time (s)');
ylabel('arithmitis theta (rad/s)');
legend('arithmitis theta', 'starithmitis theta');
grid on;

sgtitle('arithmitis and starithmitis');
%% trives

trivi = zeros(3, length(time));
strivi = zeros(3, length(time));

for i = 1:length(paronomastis)
    trivi(:, i) = inv(L)*dp(:, i)/sqrt(dp(:, i)'*inv(L)*dp(:, i));
    strivi(:, i) = inv(L)*dp_star(:, i)/sqrt(dp_star(:, i)'*inv(L)*dp_star(:, i));
end

figure;

subplot(3,1,1);
plot(time, trivi(1,:), 'LineWidth', 2);
hold on;
plot(time, strivi(1,:), 'LineWidth', 2);
title('Ground Friction x and Ground Friction Star x');
xlabel('Time (s)');
ylabel('trivi x (N)');
legend('trivi x', 'strivi x');
grid on;

subplot(3,1,2);
plot(time, trivi(2,:), 'LineWidth', 2);
hold on;
plot(time, strivi(2,:), 'LineWidth', 2);
title('Ground Friction y and Ground Friction Star y');
xlabel('Time (s)');
ylabel('trivi y (N)');
legend('trivi y', 'strivi y');
grid on;

subplot(3,1,3);
plot(time, trivi(3,:), 'LineWidth', 2);
hold on;
plot(time, strivi(3,:), 'LineWidth', 2);
title('Ground Friction theta and Ground Friction Star theta');
xlabel('Time (s)');
ylabel('trivi theta (rad/s)');
legend('strivi theta', 'strivi theta');
grid on;

sgtitle('trives anamesa se Simulated and Star');

%% Norm of velocities and friction forces

norm_trivi = zeros(size(time));
norm_strivi = zeros(size(time));
norm_dp = zeros(size(time));
norm_dp_star = zeros(size(time));

for i = 1:length(time)
    norm_trivi(i) = norm(trivi(:,i));
    norm_strivi(i) = norm(strivi(:,i));
    norm_dp(i) = norm(dp(:,i));
    norm_dp_star(i) = norm(dp_star(:,i));
end

figure;

subplot(2,1,1);
plot(time, norm_trivi, 'LineWidth', 2);
hold on;
plot(time, norm_strivi, 'LineWidth', 2);
title('Norm of friction force');
xlabel('Time (s)');
ylabel('norm of trivi (N)');
legend('norm trivi', 'norm strivi');
grid on;

subplot(2,1,2);
plot(time, norm_dp, 'LineWidth', 2);
hold on;
plot(time, norm_dp_star, 'LineWidth', 2);
title('Norm of velocities');
xlabel('Time (s)');
ylabel('Norm velocities (m/s)');
legend('norm dp', 'norm dp star');
grid on;

sgtitle('Norm of trivi and velocities');