%=======================================================================%
% This function creates a semi of the cycle trajectory,
% given:
%   - duration (Total time for the trajectory)
%   - trajectory_radius
%   - timestep
%=======================================================================%

function [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, t] = ...
                        s_shape_trajectory(duration, trajectory_radius, timestep)
    
    half_duration = duration/2;

    t1 = 0:timestep:half_duration;

    A = [1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 2 0 0 0;
         1 half_duration half_duration^2 half_duration^3 half_duration^4 half_duration^5;
         0 1 2*half_duration 3*half_duration^2 4*half_duration^3 5*half_duration^4;
         0 0 2 6*half_duration 12*half_duration^2 20*half_duration^3];

    b = [0; 0; 0; half_duration; 0; 0];
    coeffs = A \ b;
    
    [a0, a1, a2, a3, a4, a5] = deal(coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6));
    var = a0 + a1*t1 + a2*t1.^2 + a3*t1.^3 + a4*t1.^4 + a5*t1.^5;
    var_dot = a1 + 2*a2*t1 + 3*a3*t1.^2 + 4*a4*t1.^3 + 5*a5*t1.^4;
    var_ddot = 2*a2 + 6*a3*t1 + 12*a4*t1.^2 + 20*a5*t1.^3;

    x1 = trajectory_radius * cos(pi * (var/half_duration));
    y1 = trajectory_radius * sin(pi * (var/half_duration));

    x1_dot = - var_dot .* (pi/half_duration) .* trajectory_radius .* sin(pi .* (var./half_duration));
    y1_dot = var_dot .* (pi/half_duration) .* trajectory_radius .* cos(pi .* (var./half_duration));

    x1_ddot = - var_ddot .* (pi/half_duration) .* trajectory_radius .* cos(pi .* (var./half_duration));
    y1_ddot = - var_ddot .* (pi/half_duration) .* trajectory_radius .* sin(pi .* (var./half_duration));

    for i = 1:length(t1)
        theta1(i) = atan2(y1_dot(i), x1_dot(i)) - pi/2;
        theta1(i) = mod(theta1(i), 2*pi);
        if theta1(i) < 0
            theta1(i) = theta1(i) + 2*pi;
        end
    end
    theta1(1) = 0;
    theta1(length(t1)) = pi;
    % denominator = (x_star_dot.^2 + y_star_dot.^2);
    % theta_star_dot = (x_star_dot .* y_star_ddot - y_star_dot .* x_star_ddot) ./ max(denominator, eps);
    for i = 1:length(t1)-1
        theta1_dot(i) = (theta1(i+1) - theta1(i)) / timestep;
    end
    theta1_dot(length(t1)) = 0;

    t2 = half_duration:timestep:duration;

    x2 = trajectory_radius * cos(pi * (var/half_duration)) - 2*trajectory_radius;
    y2 = - trajectory_radius * sin(pi * (var/half_duration));

    x2_dot = - var_dot .* (pi/half_duration) .* trajectory_radius .* sin(pi .* (var./half_duration));
    y2_dot = - var_dot .* (pi/half_duration) .* trajectory_radius .* cos(pi .* (var./half_duration));

    x2_ddot = - var_ddot .* (pi/half_duration) .* trajectory_radius .* cos(pi .* (var./half_duration));
    y2_ddot = var_ddot .* (pi/half_duration) .* trajectory_radius .* sin(pi .* (var./half_duration));

    for i = 1:length(t2)
        theta2(i) = atan2(y2_dot(i), x2_dot(i)) - pi/2;
        theta2(i) = mod(theta2(i), 2*pi);
        if theta2(i) < 0
            theta2(i) = theta2(i) + 2*pi;
        end
    end
    theta2(1) = pi;
    theta2(length(t2)) = 0;
    % denominator = (x_star_dot.^2 + y_star_dot.^2);
    % theta_star_dot = (x_star_dot .* y_star_ddot - y_star_dot .* x_star_ddot) ./ max(denominator, eps);
    for i = 1:length(t2)-1
        theta2_dot(i) = (theta2(i+1) - theta2(i)) / timestep;
    end
    theta2_dot(length(t2)) = 0;

    x_star = [x1, x2(2:length(t2))];
    y_star = [y1, y2(2:length(t2))];
    x_star_dot = [x1_dot, x2_dot(2:length(t2))];
    y_star_dot = [y1_dot, y2_dot(2:length(t2))];
    theta_star = [theta1, theta2(2:length(t2))];
    theta_star_dot = [theta1_dot, theta2_dot(2:length(t2))];
    t = [t1, t2(2:length(t2))];
end
