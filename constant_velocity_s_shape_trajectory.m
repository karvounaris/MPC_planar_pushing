function [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, t, duration] = ...
                    constant_velocity_s_shape_trajectory(trajectory_radius, v_constant, timestep, x_center, y_center)
    % trajectory_radius: radius of the circle
    % v_constant: constant velocity
    % timestep: time step
    % x_center, y_center: center of the first semicircular trajectory
    
    % Duration
    if mod(2*pi * trajectory_radius, v_constant) == 0
        duration = floor(2*pi * trajectory_radius / v_constant);
    else
        duration = floor(2*pi * trajectory_radius / v_constant)+1;
    end
    
    % Time vector
    t = 0:timestep:duration;

    % First semicircular arc (shifted to be centered at (x_center, y_center))
    theta_star_1 = (v_constant / trajectory_radius) * t(1:floor(length(t)/2));
    theta_star_dot_1 = (v_constant / trajectory_radius) * ones(size(t(1:floor(length(t)/2)))); % Constant angular velocity
    x_star_1 = x_center + trajectory_radius * cos(theta_star_1);
    y_star_1 = y_center + trajectory_radius * sin(theta_star_1);
    x_star_dot_1 = -v_constant * sin(theta_star_1); % Derivative of x = R*cos(theta)
    y_star_dot_1 = v_constant * cos(theta_star_1);  % Derivative of y = R*sin(theta)

    % Second semicircular arc
    theta_star_2 = pi - (v_constant / trajectory_radius) * t(1:floor(length(t)/2)+1);
    theta_star_dot_2 = - (v_constant / trajectory_radius) * ones(size(t(floor(length(t)/2)+1:end))); % Constant angular velocity
    % The second arc is computed by shifting the first arc by (-2*trajectory_radius, 0)
    x_star_2 = x_center - 2*trajectory_radius + trajectory_radius * cos(pi-theta_star_2);
    y_star_2 = y_center - trajectory_radius * sin(pi-theta_star_2);
    x_star_dot_2 = -v_constant * sin(pi-theta_star_2);  % Derivative remains the same
    y_star_dot_2 = -v_constant * cos(pi-theta_star_2);  % Derivative remains the same

    % Concatenate the two arcs
    x_star = [x_star_1 x_star_2];
    y_star = [y_star_1 y_star_2];
    theta_star = [theta_star_1 theta_star_2];
    x_star_dot = [x_star_dot_1 x_star_dot_2];
    y_star_dot = [y_star_dot_1 y_star_dot_2];
    theta_star_dot = [theta_star_dot_1 theta_star_dot_2];
end
