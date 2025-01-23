function [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, phi_star, phi_star_dot, t, duration] = ...
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

    % Compute the angular acceleration (theta_ddot) as a vector
    theta_star_ddot = zeros(size(t));

    % Compute the x and y coordinates
    x_star = trajectory_radius * cos(theta_star);
    y_star = trajectory_radius * sin(theta_star);

    % Compute the velocities (x_dot and y_dot)
    x_star_dot = -v_constant * sin(theta_star); % Derivative of x = R*cos(theta)
    y_star_dot = v_constant * cos(theta_star);  % Derivative of y = R*sin(theta)

    % Compute the velocities (x_dot and y_dot)
    x_star_ddot = -trajectory_radius .* theta_star_ddot .* sin(theta_star) - trajectory_radius .* theta_star_dot.^2 .* cos(theta_star); % Derivative of x_dot = -R*theta_dot*sin(theta)
    y_star_ddot = trajectory_radius .* theta_star_ddot .* cos(theta_star) - trajectory_radius .* theta_star_dot.^2 .* sin(theta_star);  % Derivative of y_dot = R*theta_dot*cos(theta)

    % Loop over each time step to calculate phi_star
    for i = 1:length(t)
        % Extract the current values of theta_star, x_star, and y_star
        theta = theta_star(i);
        x = x_star(i);
        y = y_star(i);
        theta_dot = theta_star_dot(i);
        x_dot = x_star_dot(i);
        y_dot = y_star_dot(i);
        
        % Construct the transformation matrix gamma for the current time step
        gamma = [cos(theta) sin(theta) 0 0 0 -y*cos(theta) + x*sin(theta);
                 -sin(theta) cos(theta) 0 0 0 y*sin(theta) + x*cos(theta);
                 0 0 1 -y -x 0;
                 0 0 0 cos(theta) sin(theta) 0;
                 0 0 0 -sin(theta) cos(theta) 0;
                 0 0 0 0 0 1];

        gamma_dot = [-theta_dot*sin(theta) theta_dot*cos(theta) 0 0 0 -y_dot*cos(theta) + y*theta_dot*sin(theta) + x_dot*sin(theta) + x*theta_dot*cos(theta);
                     -theta_dot*cos(theta) -theta_dot*sin(theta) 0 0 0 y_dot*sin(theta) + y*theta_dot*cos(theta) + x_dot*cos(theta) - x*theta_dot*sin(theta);
                     0 0 0 -y_dot -x_dot 0;
                     0 0 0 -theta_dot*sin(theta) theta_dot*cos(theta) 0;
                     0 0 0 -theta_dot*cos(theta) -theta_dot*sin(theta) 0;
                     0 0 0 0 0 0];

        % Define the velocity vector
        v = [x_dot; y_dot; 0; 0; 0; theta_dot];
        
        % Define the acceleration vector
        a = [x_star_ddot(i); y_star_ddot(i); 0; 0; 0; theta_star_ddot(i)];

        % Apply the transformation to v and a
        a_tranformed = gamma_dot*v + gamma*a;
        v_tranformed = gamma*v;

        % Extract the first two components of the resulting velocity vector
        v_temp = -v_tranformed(1:2);
        a_temp = -a_tranformed(1:2);

        % Compute phi_star (angle between velocity vector and x-axis)
        phi_star(i) = atan2(v_temp(2), v_temp(1));  % Angle of the velocity vector with the x-axis

        phi_star_dot(i) = (v_temp' * [a_temp(2); -a_temp(1)]) / sqrt(v_temp(1)^2 + v_temp(2)^2);
        if phi_star_dot(i) < 1e-5
            phi_star_dot(i) = 0;
        end

        % Ensure phi_star is positive by adding 2*pi if the angle is negative
        if phi_star(i) < 0
            phi_star(i) = phi_star(i) + 2*pi;
        end
    end
end
