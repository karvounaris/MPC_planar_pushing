%=======================================================================%
% This function calculates the x and y coordinates of the contact point in
% the limit surface of the object relative to the body frame of the object,
% given:
%     - phi (angle between the contact point and the x-axis of body frame)
%     - radius
%     - len
%=======================================================================%

function [x_c, y_c, r_c, n_c, t_c] = calculate_r_c(phi, len, radius)

    phi = abs(mod(phi, 2*pi));

    if phi >= 0 && phi < pi/4
        x_c = len/2 + radius * cos(phi*2);
        y_c = radius * sin(phi*2);
        n_c = [-cos(2*phi); -sin(2*phi)];
        t_c = [sin(2*phi); -cos(2*phi)];
    elseif phi >= pi/4 && phi < 3*pi/4
        x_c = radius * (cos(phi)/sin(phi));
        y_c = radius;
        n_c = [0; -1];
        t_c = [1; 0];
    elseif phi >= 3*pi/4 && phi < 5*pi/4
        x_c = -(len/2 + radius * cos(phi*2));
        y_c = -radius * sin(phi*2);
        n_c = [cos(2*phi); sin(2*phi)];
        t_c = [-sin(2*phi); cos(2*phi)];
    elseif phi >= 5*pi/4 && phi < 7*pi/4
        x_c = -radius * (cos(phi)/sin(phi));
        y_c = -radius;
        n_c = [0; 1];
        t_c = [-1; 0];
    elseif phi >= 7*pi/4 && phi < 2*pi
        x_c = len/2 + radius * cos(phi*2);
        y_c = radius * sin(phi*2);
        n_c = [-cos(2*phi); -sin(2*phi)];
        t_c = [sin(2*phi); -cos(2*phi)];
    else
        disp("Invalid value for phi")
        x_c = 0;
        y_c = 0;
        n_c = [0; 0];
        t_c = [0; 0];
    end

    r_c = sqrt(x_c^2 + y_c^2);

end


