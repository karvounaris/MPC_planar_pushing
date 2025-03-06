function [A, B] = jacobian_function_system_linearization(x_star, u_star, L, radius, len, wid, object_shape)
    % Extract elements from x_star
    x = x_star(1);
    y = x_star(2);
    theta = x_star(3);
    phi = x_star(4);

    % Extract elements from u_star
    fn = u_star(1);
    ft = u_star(2);
    phi_dot = u_star(3);

    l1 = L(1,1);
    l2 = L(2,2);
    l3 = L(3,3);

    phi_val = x_star(4);

    if object_shape == "rectangular_capsule_prism"

        if (phi_val >= 0 && phi_val < pi/4) || (phi_val >= 7*pi/4 && phi_val < 2*pi)
            A = [0, 0, fn*(l1*cos(2*phi)*sin(theta) + l2*sin(2*phi)*cos(theta)) + ft*(l2*cos(2*phi)*cos(theta) - l1*sin(2*phi)*sin(theta)), ...
                 fn*(2*l1*sin(2*phi)*cos(theta) + 2*l2*cos(2*phi)*sin(theta)) + ft*(2*l1*cos(2*phi)*cos(theta) - 2*l2*sin(2*phi)*sin(theta));
                 0, 0, ft*(l1*sin(2*phi)*cos(theta) + l2*cos(2*phi)*sin(theta)) - fn*(l1*cos(2*phi)*cos(theta) - l2*sin(2*phi)*sin(theta)), ...
                 ft*(2*l1*cos(2*phi)*sin(theta) + 2*l2*sin(2*phi)*cos(theta)) - fn*(2*l2*cos(2*phi)*cos(theta) - 2*l1*sin(2*phi)*sin(theta));
                 0, 0, 0, ...
                 ft*l3*(2*sin(2*phi)*(len/2 + radius*cos(2*phi)) - 2*radius*cos(2*phi)*sin(2*phi)) + fn*l3*(2*radius*cos(2*phi)^2 - 2*cos(2*phi)*(len/2 + radius*cos(2*phi)));
                 0, 0, 0, 0];
            
            B = [l2*sin(2*phi)*sin(theta) - l1*cos(2*phi)*cos(theta), l1*sin(2*phi)*cos(theta) + l2*cos(2*phi)*sin(theta), 0;
                 -l1*cos(2*phi)*sin(theta) - l2*sin(2*phi)*cos(theta), l1*sin(2*phi)*sin(theta) - l2*cos(2*phi)*cos(theta), 0;
                 -l3*(sin(2*phi)*(len/2 + radius*cos(2*phi)) - radius*cos(2*phi)*sin(2*phi)), -l3*(radius*sin(2*phi)^2 + cos(2*phi)*(len/2 + radius*cos(2*phi))), 0;
                 0, 0, 1];
    
        elseif phi_val >= pi/4 && phi_val < 3*pi/4
            % Matrix A
            A = [0, 0, fn*l2*cos(theta) - ft*l1*sin(theta), 0;
                 0, 0, ft*l1*cos(theta) + fn*l2*sin(theta), 0;
                 0, 0, 0, fn*l3*radius + (fn*l3*radius*cos(phi)^2)/sin(phi)^2;
                 0, 0, 0, 0];
            
            % Matrix B
            B = [l2*sin(theta), l1*cos(theta), 0;
                 -l2*cos(theta), l1*sin(theta), 0;
                 -(l3*radius*cos(phi))/sin(phi), -l3*radius, 0;
                 0, 0, 1];
    
        elseif phi_val >= 3*pi/4 && phi_val < 5*pi/4
            % Matrix A
            A = [0, 0, - fn*(l1*cos(2*phi)*sin(theta) + l2*sin(2*phi)*cos(theta)) - ft*(l2*cos(2*phi)*cos(theta) - l1*sin(2*phi)*sin(theta)), ...
                 - fn*(2*l1*sin(2*phi)*cos(theta) + 2*l2*cos(2*phi)*sin(theta)) - ft*(2*l1*cos(2*phi)*cos(theta) - 2*l2*sin(2*phi)*sin(theta));
                 0, 0, fn*(l1*cos(2*phi)*cos(theta) - l2*sin(2*phi)*sin(theta)) - ft*(l1*sin(2*phi)*cos(theta) + l2*cos(2*phi)*sin(theta)), ...
                 fn*(2*l2*cos(2*phi)*cos(theta) - 2*l1*sin(2*phi)*sin(theta)) - ft*(2*l1*cos(2*phi)*sin(theta) + 2*l2*sin(2*phi)*cos(theta));
                 0, 0, 0, ft*l3*(2*sin(2*phi)*(len/2 + radius*cos(2*phi)) - 2*radius*cos(2*phi)*sin(2*phi)) + fn*l3*(2*radius*cos(2*phi)^2 - 2*cos(2*phi)*(len/2 + radius*cos(2*phi)));
                 0, 0, 0, 0];
            
            % Matrix B
            B = [l1*cos(2*phi)*cos(theta) - l2*sin(2*phi)*sin(theta), -l1*sin(2*phi)*cos(theta) - l2*cos(2*phi)*sin(theta), 0;
                 l1*cos(2*phi)*sin(theta) + l2*sin(2*phi)*cos(theta), l2*cos(2*phi)*cos(theta) - l1*sin(2*phi)*sin(theta), 0;
                 -l3*(sin(2*phi)*(len/2 + radius*cos(2*phi)) - radius*cos(2*phi)*sin(2*phi)), -l3*(radius*sin(2*phi)^2 + cos(2*phi)*(len/2 + radius*cos(2*phi))), 0;
                 0, 0, 1];
    
        elseif phi_val >= 5*pi/4 && phi_val < 7*pi/4
            % Matrix A
            A = [
                0, 0, ft*l1*sin(theta) - fn*l2*cos(theta), 0;
                0, 0, -ft*l1*cos(theta) - fn*l2*sin(theta), 0;
                0, 0, 0, fn*l3*radius + (fn*l3*radius*cos(phi)^2)/sin(phi)^2;
                0, 0, 0, 0];
            
            % Matrix B
            B = [-l2*sin(theta), -l1*cos(theta), 0;
                 l2*cos(theta), -l1*sin(theta), 0;
                 -(l3*radius*cos(phi))/sin(phi), -l3*radius, 0;
                 0, 0, 1];
        end

    elseif object_shape == "rectangular_prism"

        alpha = atan(wid/len);
        if ((phi_val >= 0) && (phi_val < alpha)) || ((phi_val >= 2*pi - alpha) && (phi_val < 2*pi))
            A =[0, 0, ft*l2*cos(theta) + fn*l1*sin(theta), 0;
                0, 0, ft*l2*sin(theta) - fn*l1*cos(theta), 0;
                0, 0, 0, (fn*l3*len*(tan(phi)^2 + 1))/2;
                0, 0, 0, 0];
    
            B = [-l1*cos(theta), l2*sin(theta), 0;
                 -l1*sin(theta), -l2*cos(theta), 0;
                 (l3*len*tan(phi))/2, -(l3*len)/2, 0;
                 0, 0, 1];
    
        elseif (phi_val >= alpha) && (phi_val < pi - alpha)
            % Matrix A
            A = [0, 0, fn*l2*cos(theta) - ft*l1*sin(theta), 0;
                 0, 0, ft*l1*cos(theta) + fn*l2*sin(theta), 0;
                 0, 0, 0, (fn*wid)/2 + (fn*wid*cos(phi)^2)/(2*sin(phi)^2);
                 0, 0, 0, 0];
            
            % Matrix B
            B = [l2*sin(theta), l1*cos(theta), 0;
                -l2*cos(theta), l1*sin(theta), 0;
                -(wid*cos(phi))/(2*sin(phi)), -(l3*wid)/2, 0;
                0, 0, 1];
    
        elseif (phi_val >= pi - alpha) && (phi_val < pi + alpha)
            % Matrix A
            A = [0, 0, - ft*l2*cos(theta) - fn*l1*sin(theta), 0;
                 0, 0, fn*l1*cos(theta) - ft*l2*sin(theta), 0;
                 0, 0, 0, (fn*l3*len*(tan(phi)^2 + 1))/2;
                 0, 0, 0, 0];
            
            % Matrix B
            B = [l1*cos(theta), -l2*sin(theta), 0;
                 l1*sin(theta), l2*cos(theta), 0;
                 (l3*len*tan(phi))/2, -(l3*len)/2, 0;
                 0, 0, 1];
    
        elseif (phi_val >= pi + alpha) && (phi_val < 2*pi - alpha)
            % Matrix A
            A = [0, 0, ft*l1*sin(theta) - fn*l2*cos(theta), 0;
                 0, 0, - ft*l1*cos(theta) - fn*l2*sin(theta), 0;
                 0, 0, 0, (fn*wid)/2 + (fn*wid*cos(phi)^2)/(2*sin(phi)^2);
                 0, 0, 0, 0];
            
            % Matrix B
            B = [-l2*sin(theta), -l1*cos(theta), 0;
                 l2*cos(theta), -l1*sin(theta), 0;
                 -(wid*cos(phi))/(2*sin(phi)), -(l3*wid)/2, 0;
                 0, 0, 1];
        end
    end
end