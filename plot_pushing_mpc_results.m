%% Script to import ROS2 results and replicate simulation plots

close all;
clear;
clc;

% -- 1. Read CSV File --
% data = readmatrix('/home/arl/panagiotis/panagiotis_ws/src/ROS_2_pushing_mpc_package/experiments/results/pushing_mpc_out.csv');
data = readmatrix('/home/karvounaris/Documents/diplomatiki/matlab/results/straight_line_experiment_best.csv');

time        = data(:, 1)';% 1 x N
x           = data(:, 2:5)';% 4 x N
x_star      = data(:, 6:9)';% 4 x N
x_dot       = data(:, 10:13)';% 4 x N
dx          = data(:, 14:17)';% 4 x N
u           = data(:, 18:20)';% 3 x N
u_star      = data(:, 21:23)';% 3 x N
du          = data(:, 24:26)';% 3 x N
z           = data(:, 27:29)';% 3 x N
solve_time  = data(:, 30)';% 1 x N
command_v   = data(:, 31:32)';% 2 x N
gurobi_solve_time  = data(:, 33)';% 1 x N

len = 0.2;
wid = 0.15;
radius = 0.5;
height = 0.18;
object_shape = "rectangular_prism";

%% Plot Capsule Shape Along Trajectory with Contact Points
% figure;
% plot(x(1,:), x(2,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Trajectory');
% % title('Trajectory with Capsule Shape');
% xlabel('x (m)');
% ylabel('y (m)');
% grid on;
% axis equal;
% 
% hold on;
% 
% plot(x(1,1), x(2,1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
% plot(x(1,end), x(2,end), 'yo', 'MarkerFaceColor', 'y', 'DisplayName', 'End');
% 
% if object_shape == "rectangular_capsule_prism"
%     shape_handle_handle = [];
%     % Plot the object shape at several points along the trajectory
%     for i = 1:round(length(x(1,:))/10):length(x(1,:))
%         shape_handle = get_capsule_shape(len, radius, x(1,i), x(2,i), x(3,i));
% 
%         if isempty(shape_handle_handle) 
%             shape_handle_handle = fill(shape_handle(:,1), shape_handle(:,2), 'b', 'FaceAlpha', 0.2, 'DisplayName', 'Capsule Shape');
%         else
%             fill(shape_handle(:,1), shape_handle(:,2), 'b', 'FaceAlpha', 0.2);
%         end
%     end
% 
% elseif object_shape == "rectangular_prism"
%     shape_handle_handle = [];
%     % Plot the object shape at several points along the trajectory
%     for i = 1:round(length(x(1,:))/10):length(x(1,:))
%         shape_handle = get_rectangle_shape(len, wid, x(1,i), x(2,i), x(3,i));
% 
%         if isempty(shape_handle_handle) 
%             shape_handle_handle = fill(shape_handle(:,1), shape_handle(:,2), 'b', 'FaceAlpha', 0.2, 'DisplayName', 'Capsule Shape');
%         else
%             fill(shape_handle(:,1), shape_handle(:,2), 'b', 'FaceAlpha', 0.2);
%         end
%     end
% end
% 
% contact_x_world = [];
% contact_y_world = [];
% 
% % Loop to plot contact points and unit vectors separately
% for i = 1:length(x(1,:))
%     [x_c, y_c, ~, n_c, t_c] = calculate_r_c(x(4,i), len, radius, wid, object_shape);
% 
%     R = [cos(x(3,i)), -sin(x(3,i)); sin(x(3,i)), cos(x(3,i))];
%     contact_point_global = R * [x_c; y_c] + [x(1,i); x(2,i)];
% 
%     contact_x_world = [contact_x_world; contact_point_global(1)];
%     contact_y_world = [contact_y_world; contact_point_global(2)];
% 
% end
% 
% % Plot for j == 1 (red circles)
% % plot(contact_x_world(idx_1), contact_y_world(idx_1), 'ro', 'MarkerSize', 2, 'DisplayName', 'Contact Point');
% contact_handle = plot(contact_x_world, contact_y_world, 'r-', 'LineWidth', 2, 'DisplayName', 'Contact Point');
% 
% % Plot for j == 0 (green circles)
% % plot(contact_x_world(idx_0), contact_y_world(idx_0), 'go', 'MarkerSize', 2);
% % plot(contact_x_world(idx_0), contact_y_world(idx_0), 'ro', 'MarkerSize', 2);
% 
% % plot(contact_x, contact_y, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Contact Path');
% 
% % Plot x_star and y_star trajectory
% plot(x_star(1,:), x_star(2,:), 'k-', 'LineWidth', 2, 'DisplayName', 'Desired trajectory'); % x_star trajectory in black
% plot(x_star(1,1), x_star(2,1), 'go', 'MarkerFaceColor', 'g'); % Start point of x_star in green
% plot(x_star(1,end), x_star(2,end), 'yo', 'MarkerFaceColor', 'y'); % End point of x_star in yellow
% 
% legend([shape_handle_handle; findobj(gca, 'DisplayName', 'Trajectory'); ...
%         findobj(gca, 'DisplayName', 'Desired trajectory'); ...
%         findobj(gca, 'DisplayName', 'Start'); findobj(gca, 'DisplayName', 'End'); ...
%         contact_handle], 'Location', 'Best');
% 
% hold off;

%% Plot x-x_star, y-y_star and theta-theta_star

% figure;
% 
% % First subplot for x and x_star over time
% subplot(3, 1, 1);
% plot(time, x(1, :), 'b-', 'LineWidth', 2); % Plot x
% hold on;
% plot(time, x_star(1, 1:length(x(1,:))), 'r-', 'LineWidth', 2); % Plot x_star
% % title('x and x^* Over Time');
% xlabel('Time (s)');
% ylabel('x (m)');
% legend('x', 'x^*', 'Location', 'best');
% grid on;
% 
% % Second subplot for y and y_star over time
% subplot(3, 1, 2);
% plot(time, x(2, :), 'b-', 'LineWidth', 2); % Plot y
% hold on;
% plot(time, x_star(2, 1:length(x(2,:))), 'r-', 'LineWidth', 2); % Plot y_star
% % title('y and y^* Over Time');
% xlabel('Time (s)');
% ylabel('y (m)');
% legend('y', 'y^*', 'Location', 'best');
% grid on;
% 
% % Third subplot for theta and theta_star over time
% subplot(3, 1, 3);
% plot(time, x(3, :), 'b-', 'LineWidth', 2); % Plot theta
% hold on;
% plot(time, x_star(3, 1:length(x(2,:))), 'r-', 'LineWidth', 2); % Plot theta_star
% % title('\theta and \theta^* Over Time');
% xlabel('Time (s)');
% ylabel('\theta (rad)');
% legend('\theta', '\theta^*', 'Location', 'best');
% grid on;
% 
% %% Plot dx dy and dtheta in seperate plots
% 
figure;
% First subplot for dx over time
subplot(3, 1, 1);
plot(time(1:end), dx(1,:), 'LineWidth', 2);
% title('dx Over Time');
xlabel('Time (s)');
ylabel('$\bar{x}$ (m)', 'Interpreter', 'latex');
grid on;

% Second subplot for dy over time
subplot(3, 1, 2);
plot(time(1:end), dx(2,:), 'LineWidth', 2);
% title('dy Over Time');
xlabel('Time (s)');
ylabel('$\bar{y}$ (m)', 'Interpreter', 'latex');
grid on;

% Third subplot for dtheta over time
subplot(3, 1, 3);
plot(time(1:end), dx(3,:), 'LineWidth', 2);
% title('dtheta Over Time');
xlabel('Time (s)');
ylabel('$\bar{\theta}$ (rad)', 'Interpreter', 'latex');
grid on;

% %% Plot fn ft and phi_dot in seperate plots
% 
% figure;
% % First subplot for fn over time
% subplot(3, 1, 1);
% plot(time(1:end), u(1,:), 'LineWidth', 2);
% % title('fn Over Time');
% xlabel('Time (s)');
% ylabel('fn (N)');
% grid on;
% 
% % Second subplot for ft over time
% subplot(3, 1, 2);
% plot(time(1:end), u(2,:), 'LineWidth', 2);
% % title('ft Over Time');
% xlabel('Time (s)');
% ylabel('ft (N)');
% grid on;
% 
% % Third subplot for phi_dot over time
% subplot(3, 1, 3);
% plot(time(1:end), u(3,:), 'LineWidth', 2);
% % title('phi dot Over Time');
% xlabel('Time (s)');
% ylabel('$\dot{\phi}$ (rad/sec)', 'Interpreter', 'latex');
% grid on;
% 
% 
% %% Plot dfn dft and dphi_dot in seperate plots
% 
% figure;
% 
% % First subplot for dfn over time
% subplot(3, 1, 1);
% plot(time(1:length(du(1,:))), du(1,:), 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('$\bar{f}_n$ (N)', 'Interpreter', 'latex');
% grid on;
% 
% % Second subplot for dft over time
% subplot(3, 1, 2);
% plot(time(1:length(du(2,:))), du(2,:), 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('$\bar{f}_t$ (N)', 'Interpreter', 'latex');
% grid on;
% 
% % Third subplot for dphi_dot over time
% subplot(3, 1, 3);
% plot(time(1:length(du(3,:))), du(3,:), 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('$\dot{\bar{\phi}}$ (rad/s)', 'Interpreter', 'latex');
% grid on;
% 
% 
% %% Plot z1 z2 and z3 in seperate plots
% 
% figure;
% % First subplot for z1 over time
% subplot(3, 1, 1);
% plot(time(1:end), z(1,:), 'LineWidth', 2);
% % title('z1 mode (sticking) Over Time');
% xlabel('Time (s)');
% ylabel('z1');
% grid on;
% 
% % Second subplot for z2 over time
% subplot(3, 1, 2);
% plot(time(1:end), z(2,:), 'LineWidth', 2);
% % title('z2 mode (Sliding Left) Over Time');
% xlabel('Time (s)');
% ylabel('z2');
% grid on;
% 
% % Third subplot for z3 over time
% subplot(3, 1, 3);
% plot(time(1:end), z(3,:), 'LineWidth', 2);
% % title('z3 mode (Sliding Right) Over Time');
% xlabel('Time (s)');
% ylabel('z3');
% grid on;
% 
% %% Plot phi and phi_dot in seperate plots
% 
% figure;
% 
% % First subplot for phi over time
% subplot(2, 1, 1);
% plot(time(1:end), x(4,:), 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('$\phi$ (rad)', 'Interpreter', 'latex');
% grid on;
% 
% % Second subplot for phi_dot over time
% subplot(2, 1, 2);
% plot(time(1:end), u(3,:), 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('$\dot{\phi}$ (rad/s)', 'Interpreter', 'latex');
% grid on;
% 
% %% Plots for solver times
% 
% figure;
% plot(time, solve_time, 'LineWidth', 2);
% hold on;
% plot(time, gurobi_solve_time, 'LineWidth', 2);
% % title('Solver time');
% xlabel('Simulation Time (s)');
% ylabel('Solver Time (s)');
% legend('Time By C++', 'Time By Gurobi');
% grid on;
% 
% 
% %% Robot velocity
% figure;
% subplot(2,1,1);
% plot(time(1:end), command_v(1,:), 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('$v_x$ (m/s)', 'Interpreter', 'latex');
% grid on;
% 
% subplot(2,1,2);
% plot(time(1:end), command_v(2,:), 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('$v_y$ (m/s)', 'Interpreter', 'latex');
% grid on;
% 
% %% Path error metrics

% [path_error, x_y_error, theta_error] = path_error_MIQP(x, x_star);
% 
% %%
% figure;
% subplot(2,1,1);
% plot(time(1:end), sqrt(dx(1,:).^2 + dx(2,:).^2), 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('$\sqrt{\bar{x}^2 + \bar{y}^2}$ (m)', 'Interpreter', 'latex');
% grid on;
% 
% % Second subplot for path error
% subplot(2,1,2);
% plot(time(1:end), x_y_error, 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('$\sqrt{\bar{x}^2 + \bar{y}^2}$ (m)', 'Interpreter', 'latex');
% grid on;

%% Video area
% --- Precompute contact points for each time step ---
contact_x_world = zeros(1, length(x(1,:)));
contact_y_world = zeros(1, length(x(1,:)));

for i = 1:length(x(1,:))
    [x_c, y_c, ~, ~, ~] = calculate_r_c(x(4,i), len, radius, wid, object_shape);
    R = [cos(x(3,i)), -sin(x(3,i)); sin(x(3,i)), cos(x(3,i))];
    contact_point_global = R * [x_c; y_c] + [x(1,i); x(2,i)];

    contact_x_world(i) = contact_point_global(1);
    contact_y_world(i) = contact_point_global(2);
end

% --- 1) Create and configure the video writer
videoFilename = 'experiment_best.avi';
video = VideoWriter(videoFilename);
video.FrameRate = 100;  % Adjust frame rate as desired
open(video);

% Create a figure for the animation
fig = figure();

% --- 2) Loop over each time step to create frames
for i = 1:length(x(1,:))

    % Clear figure and hold on for multiple plots
    clf;  
    hold on;  
    grid on;  
    axis equal;

    % (Optional) Set the axes to a fixed range if desired
    % xlim([-1 5]); ylim([-1 5]);  % adjust to your data

    % --- (A) Plot the completed portion of the object's trajectory so far
    plot(x(1,1:i), x(2,1:i), 'b-', 'LineWidth', 2, ...
         'DisplayName','Object Trajectory');

    % Mark the start and (current) end
    plot(x(1,1), x(2,1), 'go', 'MarkerFaceColor','g', ...
         'DisplayName','Start');
    plot(x(1,i), x(2,i), 'yo', 'MarkerFaceColor','y', ...
         'DisplayName','Current Position');

    % --- (B) Plot the desired trajectory up to the current index
    plot(x_star(1,1:i), x_star(2,1:i), 'k-', 'LineWidth', 2, ...
         'DisplayName','Desired Trajectory');

    % Mark the start and end of desired trajectory
    plot(x_star(1,1),   x_star(2,1),   'go', 'MarkerFaceColor','g');
    plot(x_star(1,end), x_star(2,end), 'yo', 'MarkerFaceColor','y');

    % --- (C) Plot the capsule (object) shape at the current step
    if object_shape == "rectangular_capsule_prism"
        object = get_capsule_shape(len, radius, x(1,i), x(2,i), x(3,i));
        fill(object(:,1), object(:,2), 'b', 'FaceAlpha', 0.2, ...
             'DisplayName','Object Shape');
    elseif object_shape == "rectangular_prism"
        object = get_rectangle_shape(len, wid, x(1,i), x(2,i), x(3,i));
        fill(object(:,1), object(:,2), 'b', 'FaceAlpha', 0.2, ...
             'DisplayName','Object Shape');
    end

    % --- (D) Plot the contact point trajectory up to current index
    % Assumes you have already computed and stored contact_x_world, contact_y_world
    % for each time step in arrays of the same length as x.
    plot(contact_x_world(1:i), contact_y_world(1:i), 'r-', 'LineWidth', 2, ...
         'DisplayName','Contact Path');

    % --- (E) Add legend and labels
    xlabel('x (m)'); ylabel('y (m)');
    % legend('Location','Best');

    % --- (F) Capture this frame and write to video
    drawnow;
    frame = getframe(fig);
    writeVideo(video, frame);
end

% --- 3) Close the video writer
close(video);
