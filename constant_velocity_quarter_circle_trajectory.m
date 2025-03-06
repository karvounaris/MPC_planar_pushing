function [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, t, duration] = ...
                    constant_velocity_quarter_circle_trajectory(trajectory_radius, v_constant, timestep, x_center, y_center)
    % trajectory_radius: radius of the circle
    % v_constant: constant velocity
    % timestep: time step
    % x_center, y_center: center of the circle

    if mod((pi/2) * trajectory_radius, v_constant) == 0
        duration = floor((pi/2) * trajectory_radius / v_constant);
    else
        duration = floor((pi/2) * trajectory_radius / v_constant)+1;
    end

    % Time vector
    t = 0:timestep:duration;

    % Compute the angle theta as a function of time
    theta_star = (v_constant / trajectory_radius) * t;

    % Compute the angular velocity (constant)
    theta_star_dot = (v_constant / trajectory_radius) * ones(size(t));

    % Compute the x and y coordinates relative to the new center
    x_star = x_center + trajectory_radius * cos(theta_star);
    y_star = y_center + trajectory_radius * sin(theta_star);

    % Compute the velocities (derivatives remain unchanged by the translation)
    x_star_dot = -v_constant * sin(theta_star);  % Derivative of x = R*cos(theta)
    y_star_dot =  v_constant * cos(theta_star);    % Derivative of y = R*sin(theta)
end