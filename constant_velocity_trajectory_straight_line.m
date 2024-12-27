%=======================================================================%
% This function creates a straight line trajectory with constant velocity,
% given:
%   - duration (Total time for the trajectory)
%   - x0, xf (Initial and final x positions)
%   - y0, yf (Initial and final y positions)
%   - theta0, thetaf (Initial and final orientations) (angles in radians)
%   - timestep
%=======================================================================%

function [x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, t] = ...
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

end
