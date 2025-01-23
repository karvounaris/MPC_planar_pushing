%=======================================================================%
% This function creates a straight line trajectory with constant velocity,
% given:
%   - duration (Total time for the trajectory)
%   - x0, xf (Initial and final x positions)
%   - y0, yf (Initial and final y positions)
%   - theta0, thetaf (Initial and final orientations) (angles in radians)
%   - timestep
%=======================================================================%

function [x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, phi_star, phi_star_dot, t] = ...
                        constant_velocity_trajectory_straight_line(duration, x_0, x_f, y_0, y_f, timestep)

    % Time vector
    t = 0:timestep:duration;
    
    % Calculate total distance and direction
    delta_x = x_f - x_0;
    delta_y = y_f - y_0;
    distance = sqrt(delta_x^2 + delta_y^2);
    
    % Calculate constant velocity components in x and y directions
    velocity = distance / duration;
    vx = velocity * (delta_x / distance); % Constant x velocity
    vy = velocity * (delta_y / distance); % Constant y velocity

    % Compute x and y positions over time
    x_star = x_0 + vx * t;
    y_star = y_0 + vy * t;

    % Velocities are constant
    x_star_dot = vx * ones(size(t));
    y_star_dot = vy * ones(size(t));
    
    % Orientation (theta) along the line
    theta_star = atan2(y_star_dot, x_star_dot) - pi/2;

    % Angular velocity (theta_dot) is zero for a straight line with constant orientation
    theta_star_dot = zeros(size(t));

    x_star_ddot = zeros(size(t));
    y_star_ddot = zeros(size(t));
    
    % Loop over each time step to calculate phi_star
    for i = 1:length(t)
        % Extract the current values of theta_star, x_star, and y_star
        theta = theta_star(i);
        x = x_star(i);
        y = y_star(i);
        
        % Construct the transformation matrix gamma for the current time step
        gamma = [cos(theta) sin(theta) 0 0 0 -y*cos(theta) + x*sin(theta);
                 -sin(theta) cos(theta) 0 0 0 y*sin(theta) + x*cos(theta);
                 0 0 1 y -x 0;
                 0 0 0 cos(theta) sin(theta) 0;
                 0 0 0 -sin(theta) cos(theta) 0;
                 0 0 0 0 0 1];

        % Define the velocity vector
        v_temp = [-x_star_dot(i); -y_star_dot(i); 0; 0; 0; theta_star_dot(i)];
        
        % Define the acceleration vector
        accel_temp = [-x_star_ddot(i); -y_star_ddot(i)];

        % Apply the transformation to v_temp
        v_temp = gamma * v_temp;

        % Extract the first two components of the resulting velocity vector
        v_temp = v_temp(1:2);

        % Compute phi_star (angle between velocity vector and x-axis)
        phi_star(i) = atan2(v_temp(2), v_temp(1));  % Angle of the velocity vector with the x-axis

        phi_star_dot(i) = (v_temp' * accel_temp) / sqrt(v_temp(1)^2 + v_temp(2)^2);

        % Ensure phi_star is positive by adding 2*pi if the angle is negative
        if phi_star(i) < 0
            phi_star(i) = phi_star(i) + 2*pi;
        end
    end
end
