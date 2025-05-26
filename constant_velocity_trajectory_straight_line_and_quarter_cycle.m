function [x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, t, duration] = ...
                        constant_velocity_trajectory_straight_line_and_quarter_cycle(v_constant, x_0, x_f, y_0, y_f, timestep, trajectory_radius)

    % Compute the difference and total distance between the start and end points
    delta_x = x_f - x_0;
    delta_y = y_f - y_0;
    distance = sqrt(delta_x^2 + delta_y^2);
    
    % Calculate the total duration based on the constant speed
    duration1 = distance / v_constant;
    duration1 = floor(duration1*1000)/1000;  % Optionally round duration to desired precision
    
    % Create the time vector
    t1 = 0:timestep:duration1;
    
    % Compute the x and y velocity components (maintaining the overall speed)
    vx = (delta_x / distance) * v_constant;
    vy = (delta_y / distance) * v_constant;
    
    % Compute positions over time ensuring the final position is reached exactly at t = duration
    x_star1 = x_0 + vx * t1;
    y_star1 = y_0 + vy * t1;
    
    % Constant velocities over time
    x_star_dot1 = vx * ones(size(t1));
    y_star_dot1 = vy * ones(size(t1));
    
    % Orientation (theta) based on the direction of travel. The subtraction of pi/2 
    % may be adjusted according to your coordinate conventions.
    theta_star1 = atan2(y_star_dot1, x_star_dot1) - pi/2;
    
    % Angular velocity is zero since the orientation remains constant along the straight line
    theta_star_dot1 = zeros(size(t1));


    x_center = x_f + trajectory_radius;
    y_center = y_f;

    if mod((pi/2) * trajectory_radius, v_constant) == 0
        duration2 = floor((pi/2) * trajectory_radius / v_constant);
    else
        duration2 = floor((pi/2) * trajectory_radius / v_constant)+1;
    end

    % Time vector
    t2 = 0:timestep:duration2;

    % Compute the angle theta as a function of time
    % theta_star2 = (v_constant / trajectory_radius) * t2 - pi/2;
    theta_star2 = -(v_constant / trajectory_radius) * t2;

    % Compute the angular velocity (constant)
    theta_star_dot2 = -(v_constant / trajectory_radius) * ones(size(t2));

    % Compute the x and y coordinates relative to the new center
    x_star2 = x_center + trajectory_radius * cos(theta_star2 - pi);
    y_star2 = y_center + trajectory_radius * sin(theta_star2 - pi);

    % Compute the velocities (derivatives remain unchanged by the translation)
    x_star_dot2 = -v_constant * sin(theta_star2);  % Derivative of x = R*cos(theta)
    y_star_dot2 =  v_constant * cos(theta_star2);    % Derivative of y = R*sin(theta)


    x_star = [x_star1 x_star2];
    y_star = [y_star1 y_star2];
    theta_star = [theta_star1 theta_star2];
    x_star_dot = [x_star_dot1 x_star_dot2];
    y_star_dot = [y_star_dot1 y_star_dot2];
    theta_star_dot = [theta_star_dot1 theta_star_dot2];
    duration = duration1 + duration2;
    t = [t1 t1(end)*ones(size(t2))+t2];
    
end

