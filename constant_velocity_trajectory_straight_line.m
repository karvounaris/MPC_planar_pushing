%=======================================================================%
% This function creates a straight line trajectory with constant velocity,
% given:
%   - duration (Total time for the trajectory)
%   - x0, xf (Initial and final x positions)
%   - y0, yf (Initial and final y positions)
%   - theta0, thetaf (Initial and final orientations) (angles in radians)
%   - timestep
%=======================================================================%

function [x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, t, duration] = ...
                        constant_velocity_trajectory_straight_line(v_constant, x_0, x_f, y_0, y_f, timestep)

    % Compute the difference and total distance between the start and end points
    delta_x = x_f - x_0;
    delta_y = y_f - y_0;
    distance = sqrt(delta_x^2 + delta_y^2);
    
    % Calculate the total duration based on the constant speed
    duration = distance / v_constant;
    duration = floor(duration*1000)/1000;  % Optionally round duration to desired precision
    
    % Create the time vector
    t = 0:timestep:duration;
    
    % Compute the x and y velocity components (maintaining the overall speed)
    vx = (delta_x / distance) * v_constant;
    vy = (delta_y / distance) * v_constant;
    
    % Compute positions over time ensuring the final position is reached exactly at t = duration
    x_star = x_0 + vx * t;
    y_star = y_0 + vy * t;
    
    % Constant velocities over time
    x_star_dot = vx * ones(size(t));
    y_star_dot = vy * ones(size(t));
    
    % Orientation (theta) based on the direction of travel. The subtraction of pi/2 
    % may be adjusted according to your coordinate conventions.
    theta_star = atan2(y_star_dot, x_star_dot) - pi/2;
    
    % Angular velocity is zero since the orientation remains constant along the straight line
    theta_star_dot = zeros(size(t));
end

