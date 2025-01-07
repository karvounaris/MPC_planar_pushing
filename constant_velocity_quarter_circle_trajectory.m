function [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, t, duration] = ...
                    constant_velocity_quarter_circle_trajectory(trajectory_radius, v_constant, timestep)
    % trajectory_radius: radius of the circle
    % v: constant velocity
    % duration: total time of the trajectory
    % timestep: time step
    
    % Duration
    duration = pi/2 * trajectory_radius / v_constant;

    % Time vector
    t = 0:timestep:duration;

    % Compute the angle theta as a function of time
    theta_star = (v_constant / trajectory_radius) * t;

    % Compute the angular velocity (theta_dot) as a vector
    theta_star_dot = (v_constant / trajectory_radius) * ones(size(t)); % Constant angular velocity

    % Compute the x and y coordinates
    x_star = trajectory_radius * cos(theta_star);
    y_star = trajectory_radius * sin(theta_star);

    % Compute the velocities (x_dot and y_dot)
    x_star_dot = -v_constant * sin(theta_star); % Derivative of x = R*cos(theta)
    y_star_dot = v_constant * cos(theta_star);  % Derivative of y = R*sin(theta)
end
