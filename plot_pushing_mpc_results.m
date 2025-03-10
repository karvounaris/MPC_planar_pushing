%% plot_pushing_mpc_results.m
% Make sure "pushing_mpc_out.csv" is in the same folder or
% specify the full path in readmatrix.

% 1) Read data from CSV
data = readmatrix('pushing_mpc_out.csv');

% Based on your CSV header:
% time, x0, x1, x2, x3, xstar0, xstar1, xstar2, xstar3,
% dx0, dx1, dx2, dx3, u0, u1, u2, du0, du1, du2, z0, z1, z2, solve_time

time       = data(:, 1);
x0         = data(:, 2);
x1         = data(:, 3);
x2         = data(:, 4);
x3         = data(:, 5);

xstar0     = data(:, 6);
xstar1     = data(:, 7);
xstar2     = data(:, 8);
xstar3     = data(:, 9);

dx0        = data(:,10);
dx1        = data(:,11);
dx2        = data(:,12);
dx3        = data(:,13);

u0         = data(:,14);
u1         = data(:,15);
u2         = data(:,16);

du0        = data(:,17);
du1        = data(:,18);
du2        = data(:,19);

z0         = data(:,20);
z1         = data(:,21);
z2         = data(:,22);

solve_time = data(:,23);

% For easier plotting in the style of your script,
% create "row-vectors" so that x(i,:) has all time samples in columns.
x      = [x0';   x1';   x2';   x3'];   % system states
x_star = [xstar0';xstar1';xstar2';xstar3']; % desired states
dx     = [dx0';  dx1';  dx2';  dx3'];  % velocities
u      = [u0';   u1';   u2'];         % control input
du     = [du0';  du1';  du2'];        % MIQP output
z      = [z0';   z1';   z2'];         % binary variables

%% ------------------------------------------------------------
%  2) PLOT: x-y Trajectory vs Desired
% -------------------------------------------------------------
figure;
plot(x(1,:), x(2,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Trajectory'); hold on;
plot(x_star(1,:), x_star(2,:), 'k--', 'LineWidth', 2, 'DisplayName', 'Desired Traj');
plot(x(1,1), x(2,1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(x(1,end), x(2,end), 'yo', 'MarkerFaceColor', 'y', 'DisplayName', 'End');
xlabel('x (m)'); ylabel('y (m)');
axis equal; grid on;
legend('Location','best');
title('x-y Trajectory vs. Desired');

%% ------------------------------------------------------------
%  3) PLOT: State vs Desired over time
%     (x, y, theta, phi) vs (x*, y*, theta*, phi*)
% -------------------------------------------------------------
figure;

% x and x*
subplot(4,1,1);
plot(time, x(1,:), 'b-', 'LineWidth', 2); hold on;
plot(time, x_star(1,:), 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('x (m)');
legend('x','x^*','Location','best'); grid on;

% y and y*
subplot(4,1,2);
plot(time, x(2,:), 'b-', 'LineWidth', 2); hold on;
plot(time, x_star(2,:), 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('y (m)');
legend('y','y^*','Location','best'); grid on;

% theta and theta*
subplot(4,1,3);
plot(time, x(3,:), 'b-', 'LineWidth', 2); hold on;
plot(time, x_star(3,:), 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\theta (rad)');
legend('\theta','\theta^*','Location','best'); grid on;

% phi and phi*
subplot(4,1,4);
plot(time, x(4,:), 'b-', 'LineWidth', 2); hold on;
plot(time, x_star(4,:), 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\phi (rad)');
legend('\phi','\phi^*','Location','best'); grid on;

sgtitle('States vs Desired States');

%% ------------------------------------------------------------
%  4) PLOT: dx, dy, dtheta, dphi (the system velocities)
% -------------------------------------------------------------
figure;
subplot(4,1,1);
plot(time, dx(1,:), 'LineWidth', 2);
xlabel('Time (s)'); ylabel('$\dot{x}$ (m/s)','Interpreter','latex');
grid on;

subplot(4,1,2);
plot(time, dx(2,:), 'LineWidth', 2);
xlabel('Time (s)'); ylabel('$\dot{y}$ (m/s)','Interpreter','latex');
grid on;

subplot(4,1,3);
plot(time, dx(3,:), 'LineWidth', 2);
xlabel('Time (s)'); ylabel('$\dot{\theta}$ (rad/s)','Interpreter','latex');
grid on;

subplot(4,1,4);
plot(time, dx(4,:), 'LineWidth', 2);
xlabel('Time (s)'); ylabel('$\dot{\phi}$ (rad/s)','Interpreter','latex');
grid on;

sgtitle('System Velocities');

%% ------------------------------------------------------------
%  5) PLOT: Control inputs u0, u1, u2 over time
% -------------------------------------------------------------
figure;
subplot(3,1,1);
plot(time, u(1,:), 'LineWidth', 2);
xlabel('Time (s)'); ylabel('u0');
grid on;

subplot(3,1,2);
plot(time, u(2,:), 'LineWidth', 2);
xlabel('Time (s)'); ylabel('u1');
grid on;

subplot(3,1,3);
plot(time, u(3,:), 'LineWidth', 2);
xlabel('Time (s)'); ylabel('u2');
grid on;

sgtitle('Control Inputs (u0, u1, u2)');

%% ------------------------------------------------------------
%  6) PLOT: MIQP increments du0, du1, du2
% -------------------------------------------------------------
figure;
subplot(3,1,1);
plot(time, du(1,:), 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta u0');
grid on;

subplot(3,1,2);
plot(time, du(2,:), 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta u1');
grid on;

subplot(3,1,3);
plot(time, du(3,:), 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta u2');
grid on;

sgtitle('MIQP Updates (\Delta u)');

%% ------------------------------------------------------------
%  7) PLOT: Binary variables z0, z1, z2
% -------------------------------------------------------------
figure;
subplot(3,1,1);
plot(time, z(1,:), 'LineWidth', 2);
xlabel('Time (s)'); ylabel('z0'); grid on;

subplot(3,1,2);
plot(time, z(2,:), 'LineWidth', 2);
xlabel('Time (s)'); ylabel('z1'); grid on;

subplot(3,1,3);
plot(time, z(3,:), 'LineWidth', 2);
xlabel('Time (s)'); ylabel('z2'); grid on;

sgtitle('Binary Variables (z0, z1, z2)');

%% ------------------------------------------------------------
%  8) PLOT: Solver time
% -------------------------------------------------------------
figure;
plot(time, solve_time, 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Solve time (s)');
title('MIQP Solve Time per Timestep');
grid on;
