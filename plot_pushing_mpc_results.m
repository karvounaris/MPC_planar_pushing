%% Script to import ROS2 results and replicate simulation plots

% -- 1. Read CSV File --
data = readmatrix('pushing_mpc_out.csv');

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
height = 0.18;
object_shape = "rectangular_prism";

%% Plot 1: Trajectory + Capsule Shape
figure;
plot(x(1,:), x(2,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Trajectory');
xlabel('x (m)');
ylabel('y (m)');
grid on;
axis equal;
hold on;

% Start and end
plot(x(1,1), x(2,1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(x(1,end), x(2,end), 'yo', 'MarkerFaceColor', 'y', 'DisplayName', 'End');

% -- Rectangular shape along the path --
rect_shape_handle = [];
step = round(length(x(1,:))/10);
for i = 1:step:length(x(1,:))
    rect_vertices = get_rectangle_shape(len, wid, x(1,i), x(2,i), x(3,i));
    if isempty(rect_shape_handle)
        rect_shape_handle = fill(rect_vertices(:,1), rect_vertices(:,2), 'b', ...
                                 'FaceAlpha', 0.2, 'DisplayName', 'Rectangle Shape');
    else
        fill(rect_vertices(:,1), rect_vertices(:,2), 'b', 'FaceAlpha', 0.2);
    end
end

contact_x_world = [];
contact_y_world = [];
for i = 1:length(x(1,:))
    [x_c, y_c, ~, n_c, t_c] = calculate_r_c(x(4,i), len, radius, wid, object_shape);
    R = [cos(x(3,i)), -sin(x(3,i));
         sin(x(3,i)),  cos(x(3,i))];
    contact_point_global = R * [x_c; y_c] + [x(1,i); x(2,i)];
    contact_x_world = [contact_x_world; contact_point_global(1)];
    contact_y_world = [contact_y_world; contact_point_global(2)];
end

contact_handle = plot(contact_x_world, contact_y_world, 'r-', 'LineWidth', 2, 'DisplayName', 'Contact Point');

% Desired trajectory in black
plot(x_star(1,:), x_star(2,:), 'k-', 'LineWidth', 2, 'DisplayName', 'Desired trajectory');
plot(x_star(1,1), x_star(2,1), 'go', 'MarkerFaceColor', 'g');
plot(x_star(1,end), x_star(2,end), 'yo', 'MarkerFaceColor', 'y');

% Adjust legend as needed
legend([rect_shape_handle; findobj(gca, 'DisplayName', 'Trajectory'); ...
        findobj(gca, 'DisplayName', 'Desired trajectory'); ...
        findobj(gca, 'DisplayName', 'Start'); findobj(gca, 'DisplayName', 'End'); ...
        contact_handle], 'Location', 'Best');
hold off;

%% Plot 2: x-x*, y-y*, theta-theta* vs Time
figure;

% x / x_star
subplot(4,1,1);
plot(time, x(1,:), 'b-', 'LineWidth', 2);
hold on;
plot(time, x_star(1, :), 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('x (m)');
legend('x', 'x^*', 'Location', 'best');
grid on;

% y / y_star
subplot(4,1,2);
plot(time, x(2,:), 'b-', 'LineWidth', 2);
hold on;
plot(time, x_star(2, :), 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('y (m)');
legend('y', 'y^*', 'Location', 'best');
grid on;

% theta / theta_star
subplot(4,1,3);
plot(time, x(3,:), 'b-', 'LineWidth', 2);
hold on;
plot(time, x_star(3, :), 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('\theta (rad)');
legend('\theta', '\theta^*', 'Location', 'best');
grid on;

subplot(4,1,4);
plot(time, x(4,:), 'b-', 'LineWidth', 2);
hold on;
plot(time, x_star(4, :), 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\phi$ (rad)', 'Interpreter', 'latex');
grid on;

%% Plot 3: xdot-xdot*, ydot-ydot*, thetadot-thetadot* vs Time
figure;

% x / x_star
subplot(4,1,1);
plot(time, x_dot(1,:), 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{x}$ (m/s)', 'Interpreter', 'latex');
grid on;

% y / y_star
subplot(4,1,2);
plot(time, x_dot(2,:), 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{y}$ (m/s)', 'Interpreter', 'latex');
grid on;

% theta / theta_star
subplot(4,1,3);
plot(time, x_dot(3,:), 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{\theta}$ (m/s)', 'Interpreter', 'latex');
grid on;

subplot(4,1,4);
plot(time, x_dot(4,:), 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{\phi}$ (m/s)', 'Interpreter', 'latex');
grid on;

%% Plot 4: dx, dy, dtheta over time

figure;
subplot(3,1,1);
plot(time, dx(1,:), 'LineWidth', 2);
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('$\bar{x}$ (m)', 'Interpreter', 'latex');
grid on;

subplot(3,1,2);
plot(time, dx(2,:), 'LineWidth', 2);
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('$\bar{y}$ (m)', 'Interpreter', 'latex');
grid on;

subplot(3,1,3);
plot(time, dx(3,:), 'LineWidth', 2);
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('$\bar{\theta}$ (rad)', 'Interpreter', 'latex');
grid on;

figure;
subplot(2,1,1);
plot(time, sqrt(dx(1,:).^2 + dx(2,:).^2), 'LineWidth', 2);
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('$\sqrt{\bar{x}^2 + \bar{y}^2}$ (m)', 'Interpreter', 'latex');
grid on;

subplot(2,1,2);
plot(time, sqrt(dx(1,:).^2 + dx(2,:).^2 + dx(3,:).^2), 'LineWidth', 2);
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('$\sqrt{\bar{x}^2 + \bar{y}^2 + \bar{\theta}^2}$', 'Interpreter', 'latex');
grid on;

%% Plot 5: u = [fn, ft, phi_dot]
figure;
subplot(3,1,1);
plot(time(1:end), u(1,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$f_n$ (N)', 'Interpreter', 'latex');
grid on;

subplot(3,1,2);
plot(time(1:end), u(2,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$f_t$ (N)', 'Interpreter', 'latex');
grid on;

subplot(3,1,3);
plot(time(1:end), u(3,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{\phi}$ (rad/s)', 'Interpreter', 'latex');
grid on;


%% Plot 6: du = [dfn, dft, dphi_dot]
figure;
subplot(3,1,1);
plot(time(1:end), du(1,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\bar{f}_n$ (N)', 'Interpreter', 'latex');
grid on;

subplot(3,1,2);
plot(time(1:end), du(2,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\bar{f}_t$ (N)', 'Interpreter', 'latex');
grid on;

subplot(3,1,3);
plot(time(1:end), du(3,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\dot{\bar{\phi}}$ (rad/s)', 'Interpreter', 'latex');
grid on;


%% Plot 7: z = [z1, z2, z3]
figure;
subplot(3,1,1);
plot(time(1:end), z(1,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('z1');
grid on;

subplot(3,1,2);
plot(time(1:end), z(2,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('z2');
grid on;

subplot(3,1,3);
plot(time(1:end), z(3,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('z3');
grid on;

%% Plot 9: Solver Times
figure;
plot(time, solve_time, 'LineWidth', 2);
hold on;
plot(time, gurobi_solve_time, 'LineWidth', 2);
xlabel('Simulation Time (s)');
ylabel('Solver Time (s)');
legend('Time By C++', 'Time By Gurobi');
grid on;

%% Plot 10: robot velocity
figure;
subplot(2,1,1);
plot(time(1:end), command_v(1,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$v_x$ (m/s)', 'Interpreter', 'latex');
grid on;

subplot(2,1,2);
plot(time(1:end), command_v(2,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$v_y$ (m/s)', 'Interpreter', 'latex');
grid on;

%% Plot 11: Path Error

% [path_error, x_y_error, theta_error] = path_error_MIQP(x, x_star);
% 
% figure;
% subplot(3,1,1);
% plot(time, path_error, 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('$\sqrt{\bar{x}^2 + \bar{y}^2 + \bar{\theta}^2}$', 'Interpreter', 'latex');
% grid on;
% 
% subplot(3,1,2);
% plot(time, x_y_error, 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('$\sqrt{\bar{x}^2 + \bar{y}^2}$', 'Interpreter', 'latex');
% grid on;
% 
% subplot(3,1,3);
% plot(time, theta_error, 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('$\bar{\theta}$ (rad)', 'Interpreter', 'latex');
% grid on;
