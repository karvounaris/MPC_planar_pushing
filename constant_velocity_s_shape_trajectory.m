function [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, t, duration] = ...
                    constant_velocity_s_shape_trajectory(trajectory_radius, v_constant, timestep)
    % trajectory_radius: radius of the circle
    % v_constant: constant velocity
    % timestep: time step
    
    % Duration
    duration = 2*pi * trajectory_radius / v_constant;

    % Time vector
    t = 0:timestep:duration;

    % Compute the angle theta as a function of time
    theta_star_1 = (v_constant / trajectory_radius) * t(1:end/2);

    % Compute the angular velocity (theta_dot) as a vector
    theta_star_dot_1 = (v_constant / trajectory_radius) * ones(size(t(1:end/2))); % Constant angular velocity

    % Compute the x and y coordinates
    x_star_1 = trajectory_radius * cos(theta_star_1);
    y_star_1 = trajectory_radius * sin(theta_star_1);

    % Compute the velocities (x_dot and y_dot)
    x_star_dot_1 = -v_constant * sin(theta_star_1); % Derivative of x = R*cos(theta)
    y_star_dot_1 = v_constant * cos(theta_star_1);  % Derivative of y = R*sin(theta)

    % Compute the angle theta as a function of time
    theta_star_2 = pi - (v_constant / trajectory_radius) * t(1:end/2);

    % Compute the angular velocity (theta_dot) as a vector
    theta_star_dot_2 = -(v_constant / trajectory_radius) * ones(size(t(end/2+1:end))); % Constant angular velocity

    % Compute the x and y coordinates
    x_star_2 = -2*trajectory_radius + trajectory_radius * cos(theta_star_1);
    y_star_2 = -trajectory_radius * sin(theta_star_1);

    % Compute the velocities (x_dot and y_dot)
    x_star_dot_2 = -v_constant * sin(theta_star_1);  % Derivative of x = R*cos(theta)
    y_star_dot_2 = -v_constant * cos(theta_star_1);  % Derivative of y = R*sin(theta)

    x_star = [x_star_1 x_star_2];
    y_star = [y_star_1 y_star_2];
    theta_star = [theta_star_1 theta_star_2];
    x_star_dot = [x_star_dot_1 x_star_dot_2];
    y_star_dot = [y_star_dot_1 y_star_dot_2];
    theta_star_dot = [theta_star_dot_1 theta_star_dot_2];
end
