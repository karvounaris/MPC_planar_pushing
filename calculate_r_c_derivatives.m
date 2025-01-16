%=======================================================================%
% This function calculates the x_c and y_c partial derivatives in respect
% phi relative to the body frame of the object,
% given:
%     - phi (angle between the contact point and the x-axis of body frame)
%     - radius
%=======================================================================%

function r_c_partial_deriative_phi = calculate_r_c_derivatives(phi, radius)

    phi = abs(mod(phi, 2*pi));

    if phi >= 0 && phi < pi/4
        x_c_partial_deriative_phi = - 2 * radius * sin(2*phi);
        y_c_partial_deriative_phi = 2 * radius * cos(2*phi);
    elseif phi >= pi/4 && phi < 3*pi/4
        x_c_partial_deriative_phi = - radius * ((cos(phi)^2 / sin(phi)^2) + 1);
        y_c_partial_deriative_phi = 0;
    elseif phi >= 3*pi/4 && phi < 5*pi/4
        x_c_partial_deriative_phi = 2 * radius * sin(2*phi);
        y_c_partial_deriative_phi = - 2 * radius * cos(2*phi);
    elseif phi >= 5*pi/4 && phi < 7*pi/4
        x_c_partial_deriative_phi = radius * ((cos(phi)^2 / sin(phi)^2) + 1);
        y_c_partial_deriative_phi = 0;
    elseif phi >= 7*pi/4 && phi < 2*pi
        x_c_partial_deriative_phi = - 2 * radius * sin(2*phi);
        y_c_partial_deriative_phi = 2 * radius * cos(2*phi);        
    else
        disp("Invalid value for phi")
        x_c_partial_deriative_phi = 0;
        y_c_partial_deriative_phi = 0;
    end

    r_c_partial_deriative_phi = [x_c_partial_deriative_phi; y_c_partial_deriative_phi];

end


