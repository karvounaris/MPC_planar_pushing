function [A, B] = jacobian_function_system_linearization(x_star, u_star, L, radius, len)
    % Extract elements from x_star
    x_sym = x_star(1);
    y_sym = x_star(2);
    theta_sym = x_star(3);
    phi_sym = x_star(4);

    % Extract elements from u_star
    fn_sym = u_star(1);
    ft_sym = u_star(2);
    phi_dot_sym = u_star(3);

    l1 = L(1,1);
    l2 = L(2,2);
    l3 = L(3,3);

    phi_val = x_star(4);

    if (phi_val >= 0 && phi_val < pi/4) || (phi_val >= 7*pi/4 && phi_val < 2*pi)
        A = [0, 0, fn_sym*(l1*cos(2*phi_sym)*sin(theta_sym) + l2*sin(2*phi_sym)*cos(theta_sym)) + ft_sym*(l2*cos(2*phi_sym)*cos(theta_sym) - l1*sin(2*phi_sym)*sin(theta_sym)), ...
             fn_sym*(2*l1*sin(2*phi_sym)*cos(theta_sym) + 2*l2*cos(2*phi_sym)*sin(theta_sym)) + ft_sym*(2*l1*cos(2*phi_sym)*cos(theta_sym) - 2*l2*sin(2*phi_sym)*sin(theta_sym));
             0, 0, ft_sym*(l1*sin(2*phi_sym)*cos(theta_sym) + l2*cos(2*phi_sym)*sin(theta_sym)) - fn_sym*(l1*cos(2*phi_sym)*cos(theta_sym) - l2*sin(2*phi_sym)*sin(theta_sym)), ...
             ft_sym*(2*l1*cos(2*phi_sym)*sin(theta_sym) + 2*l2*sin(2*phi_sym)*cos(theta_sym)) - fn_sym*(2*l2*cos(2*phi_sym)*cos(theta_sym) - 2*l1*sin(2*phi_sym)*sin(theta_sym));
             0, 0, 0, ...
             ft_sym*l3*(2*sin(2*phi_sym)*(len/2 + radius*cos(2*phi_sym)) - 2*radius*cos(2*phi_sym)*sin(2*phi_sym)) + fn_sym*l3*(2*radius*cos(2*phi_sym)^2 - 2*cos(2*phi_sym)*(len/2 + radius*cos(2*phi_sym)));
             0, 0, 0, 0];
        
        B = [l2*sin(2*phi_sym)*sin(theta_sym) - l1*cos(2*phi_sym)*cos(theta_sym), l1*sin(2*phi_sym)*cos(theta_sym) + l2*cos(2*phi_sym)*sin(theta_sym), 0;
             -l1*cos(2*phi_sym)*sin(theta_sym) - l2*sin(2*phi_sym)*cos(theta_sym), l1*sin(2*phi_sym)*sin(theta_sym) - l2*cos(2*phi_sym)*cos(theta_sym), 0;
             -l3*(sin(2*phi_sym)*(len/2 + radius*cos(2*phi_sym)) - radius*cos(2*phi_sym)*sin(2*phi_sym)), -l3*(radius*sin(2*phi_sym)^2 + cos(2*phi_sym)*(len/2 + radius*cos(2*phi_sym))), 0;
             0, 0, 1];

    elseif phi_val >= pi/4 && phi_val < 3*pi/4
        % Matrix A
        A = [0, 0, fn_sym*l2*cos(theta_sym) - ft_sym*l1*sin(theta_sym), 0;
             0, 0, ft_sym*l1*cos(theta_sym) + fn_sym*l2*sin(theta_sym), 0;
             0, 0, 0, fn_sym*l3*radius + (fn_sym*l3*radius*cos(phi_sym)^2)/sin(phi_sym)^2;
             0, 0, 0, 0];
        
        % Matrix B
        B = [l2*sin(theta_sym), l1*cos(theta_sym), 0;
             -l2*cos(theta_sym), l1*sin(theta_sym), 0;
             -(l3*radius*cos(phi_sym))/sin(phi_sym), -l3*radius, 0;
             0, 0, 1];

    elseif phi_val >= 3*pi/4 && phi_val < 5*pi/4
        % Matrix A
        A = [0, 0, - fn_sym*(l1*cos(2*phi_sym)*sin(theta_sym) + l2*sin(2*phi_sym)*cos(theta_sym)) - ft_sym*(l2*cos(2*phi_sym)*cos(theta_sym) - l1*sin(2*phi_sym)*sin(theta_sym)), ...
             - fn_sym*(2*l1*sin(2*phi_sym)*cos(theta_sym) + 2*l2*cos(2*phi_sym)*sin(theta_sym)) - ft_sym*(2*l1*cos(2*phi_sym)*cos(theta_sym) - 2*l2*sin(2*phi_sym)*sin(theta_sym));
             0, 0, fn_sym*(l1*cos(2*phi_sym)*cos(theta_sym) - l2*sin(2*phi_sym)*sin(theta_sym)) - ft_sym*(l1*sin(2*phi_sym)*cos(theta_sym) + l2*cos(2*phi_sym)*sin(theta_sym)), ...
             fn_sym*(2*l2*cos(2*phi_sym)*cos(theta_sym) - 2*l1*sin(2*phi_sym)*sin(theta_sym)) - ft_sym*(2*l1*cos(2*phi_sym)*sin(theta_sym) + 2*l2*sin(2*phi_sym)*cos(theta_sym));
             0, 0, 0, ft_sym*l3*(2*sin(2*phi_sym)*(len/2 + radius*cos(2*phi_sym)) - 2*radius*cos(2*phi_sym)*sin(2*phi_sym)) + fn_sym*l3*(2*radius*cos(2*phi_sym)^2 - 2*cos(2*phi_sym)*(len/2 + radius*cos(2*phi_sym)));
             0, 0, 0, 0];
        
        % Matrix B
        B = [l1*cos(2*phi_sym)*cos(theta_sym) - l2*sin(2*phi_sym)*sin(theta_sym), -l1*sin(2*phi_sym)*cos(theta_sym) - l2*cos(2*phi_sym)*sin(theta_sym), 0;
             l1*cos(2*phi_sym)*sin(theta_sym) + l2*sin(2*phi_sym)*cos(theta_sym), l2*cos(2*phi_sym)*cos(theta_sym) - l1*sin(2*phi_sym)*sin(theta_sym), 0;
             -l3*(sin(2*phi_sym)*(len/2 + radius*cos(2*phi_sym)) - radius*cos(2*phi_sym)*sin(2*phi_sym)), -l3*(radius*sin(2*phi_sym)^2 + cos(2*phi_sym)*(len/2 + radius*cos(2*phi_sym))), 0;
             0, 0, 1];

    elseif phi_val >= 5*pi/4 && phi_val < 7*pi/4
        % Matrix A
        A = [
            0, 0, ft_sym*l1*sin(theta_sym) - fn_sym*l2*cos(theta_sym), 0;
            0, 0, -ft_sym*l1*cos(theta_sym) - fn_sym*l2*sin(theta_sym), 0;
            0, 0, 0, fn_sym*l3*radius + (fn_sym*l3*radius*cos(phi_sym)^2)/sin(phi_sym)^2;
            0, 0, 0, 0];
        
        % Matrix B
        B = [-l2*sin(theta_sym), -l1*cos(theta_sym), 0;
             l2*cos(theta_sym), -l1*sin(theta_sym), 0;
             -(l3*radius*cos(phi_sym))/sin(phi_sym), -l3*radius, 0;
             0, 0, 1];
    end
end