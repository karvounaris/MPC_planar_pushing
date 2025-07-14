function [x_c, y_c, r_c, n_c, t_c] = calculate_r_c(phi, len, radius, wid, object_shape)

    phi = abs(mod(phi, 2*pi));
    if phi == 2*pi
        phi = 0;
    end
    if object_shape == "rectangular_capsule_prism"
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

    elseif object_shape == "rectangular_prism"
        alpha = atan(wid/len);   % corner angle

        if ((phi >= 0) && (phi < alpha)) || ((phi >= 2*pi - alpha) && (phi < 2*pi))
            % Right edge
            x_c = len/2;
            y_c = (len/2)*tan(phi);
            n_c = [-1; 0];
            t_c = [0; -1];
    
        elseif (phi >= alpha) && (phi < pi - alpha)
            % Top edge
            x_c = wid/2 * (cos(phi)./sin(phi));
            y_c = wid/2;
            n_c = [0; -1];
            t_c = [1; 0];
    
        elseif (phi >= pi - alpha) && (phi < pi + alpha)
            x_c = -len/2;
            y_c = -len/2 * tan(phi);
            n_c = [1; 0];
            t_c = [0; 1];
    
        elseif (phi >= pi + alpha) && (phi < 2*pi - alpha)
            x_c = -(wid/2)*cos(phi)/sin(phi);
            y_c = -wid/2;
            n_c = [0; 1];
            t_c = [-1; 0];
    
        else
            x_c = 0; 
            y_c = 0;
            n_c = [0; 0];
            t_c = [0; 0];
        end
    
        % Once x_c, y_c are set, compute distance
        r_c = sqrt(x_c^2 + y_c^2);
    end

end


