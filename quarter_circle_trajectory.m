%=======================================================================%
% This function creates a semi of the cycle trajectory,
% given:
%   - duration (Total time for the trajectory)
%   - trajectory_radius
%   - timestep
%=======================================================================%

function [x_star, y_star, x_star_dot, y_star_dot, theta_star, theta_star_dot, t] = ...
                        quarter_circle_trajectory(duration, trajectory_radius, timestep)
    t = 0:timestep:duration;

    A = [1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 2 0 0 0;
         1 duration duration^2 duration^3 duration^4 duration^5;
         0 1 2*duration 3*duration^2 4*duration^3 5*duration^4;
         0 0 2 6*duration 12*duration^2 20*duration^3];

    b = [0; 0; 0; duration; 0; 0];
    coeffs = A \ b;

    [a0, a1, a2, a3, a4, a5] = deal(coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6));
    var = a0 + a1*t + a2*t.^2 + a3*t.^3 + a4*t.^4 + a5*t.^5;
    var_dot = a1 + 2*a2*t + 3*a3*t.^2 + 4*a4*t.^3 + 5*a5*t.^4;
    var_ddot = 2*a2 + 6*a3*t + 12*a4*t.^2 + 20*a5*t.^3;

    x_star = trajectory_radius * cos(pi/2 * (var/duration));
    y_star = trajectory_radius * sin(pi/2 * (var/duration));

    x_star_dot = - var_dot .* (pi/(2*duration)) .* trajectory_radius .* sin(pi/2 .* (var./duration));
    y_star_dot = var_dot .* (pi/(2*duration)) .* trajectory_radius .* cos(pi/2 .* (var./duration));

    x_star_ddot = - var_ddot .* (pi/(2*duration)) .* trajectory_radius .* cos(pi/2 .* (var./duration));
    y_star_ddot = - var_ddot .* (pi/(2*duration)) .* trajectory_radius .* sin(pi/2 .* (var./duration));

    for i = 1:length(t)
        theta_star(i) = atan2(y_star_dot(i), x_star_dot(i)) - pi/2;
        theta_star(i) = mod(theta_star(i), 2*pi);
        if theta_star(i) < 0
            theta_star(i) = theta_star(i) + 2*pi;
        end
    end
    theta_star(1) = 0;
    theta_star(length(t)) = pi/2;
    for i = 1:length(t)-1
        theta_star_dot(i) = (theta_star(i+1) - theta_star(i)) / timestep;
    end
    theta_star_dot(length(t)) = 0;
end




