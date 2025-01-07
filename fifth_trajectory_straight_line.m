%=======================================================================%
% This function creates a 5th-degree polynomial trajectory along the +y axis,
% given:
%   - T (Total time for the trajectory)
%   - x0, xf (Initial and final x positions)
%   - y0, yf (Initial and final y positions)
%   - theta0, thetaf (Initial and final orientations) (angles in radians)
%   - timestep
%=======================================================================%

function [x_star, x_star_dot, y_star, y_star_dot, theta_star, theta_star_dot, t] = ...
                        fifth_trajectory_straight_line(duration, x_0, x_f, y_0, y_f, timestep)

    t = 0:timestep:duration;

    A = [1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 2 0 0 0;
         1 duration duration^2 duration^3 duration^4 duration^5;
         0 1 2*duration 3*duration^2 4*duration^3 5*duration^4;
         0 0 2 6*duration 12*duration^2 20*duration^3];

    bx = [x_0; 0; 0; x_f; 0; 0];
    coeffs_x = A \ bx;

    by = [y_0; 0; 0; y_f; 0; 0];
    coeffs_y = A \ by;

    [a0x, a1x, a2x, a3x, a4x, a5x] = deal(coeffs_x(1), coeffs_x(2), coeffs_x(3), coeffs_x(4), coeffs_x(5), coeffs_x(6));
    [a0y, a1y, a2y, a3y, a4y, a5y] = deal(coeffs_y(1), coeffs_y(2), coeffs_y(3), coeffs_y(4), coeffs_y(5), coeffs_y(6));

    x_star = a0x + a1x*t + a2x*t.^2 + a3x*t.^3 + a4x*t.^4 + a5x*t.^5;
    y_star = a0y + a1y*t + a2y*t.^2 + a3y*t.^3 + a4y*t.^4 + a5y*t.^5;

    x_star_dot = a1x + 2*a2x*t + 3*a3x*t.^2 + 4*a4x*t.^3 + 5*a5x*t.^4;
    y_star_dot = a1y + 2*a2y*t + 3*a3y*t.^2 + 4*a4y*t.^3 + 5*a5y*t.^4;

    x_star_ddot = 2*a2x + 6*a3x*t + 12*a4x*t.^2 + 20*a5x*t.^3;
    y_star_ddot = 2*a2y + 6*a3y*t + 12*a4y*t.^2 + 20*a5y*t.^3;

    theta_star = atan2(y_star_dot, x_star_dot) - pi/2;
    % theta_star(1) = theta_star(1) + pi;
    theta_star(1) = theta_star(2);
    theta_star(end) = theta_star(2);

    for i = 1:length(t)-1
        theta_star_dot(i) = (theta_star(i+1) - theta_star(i)) / timestep;
    end
    theta_star_dot(length(t)) = 0;

end
