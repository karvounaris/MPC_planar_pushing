function [A, B] = jacobian_function_system_linearization_square(x_star, u_star, L, len)
    % Extract elements from x_star
    x = x_star(1);
    y = x_star(2);
    theta_sym = x_star(3);
    phi_sym = x_star(4);

    % Extract elements from u_star
    fn = u_star(1);
    ft = u_star(2);
    phi_dot = u_star(3);

    l1 = L(1,1);
    l2 = L(2,2);
    l3 = L(3,3);

    phi_val = x_star(4);

    if (phi_val >= 0 && phi_val < pi/4) || (phi_val >= 7*pi/4 && phi_val < 2*pi)
        A =[0, 0, ft*l2*cos(theta) + fn*l1*sin(theta), 0;
            0, 0, ft*l2*sin(theta) - fn*l1*cos(theta), 0;
            0, 0, 0, (fn*l3*len)/2 + (fn*l3*len*sin(phi)^2)/(2*cos(phi)^2);
            0, 0, 0, 0];

        B = [-l1*cos(theta), l2*sin(theta), 0;
             -l1*sin(theta), -l2*cos(theta), 0;
             (l3*len*sin(phi))/(2*cos(phi)), -(l3*len)/2, 0;
             0, 0, 1];

    elseif phi_val >= pi/4 && phi_val < 3*pi/4
        % Matrix A
        A = [0, 0, fn*l2*cos(theta) - ft*l1*sin(theta), 0;
             0, 0, ft*l1*cos(theta) + fn*l2*sin(theta), 0;
             0, 0, 0, (fn*l3*len)/2 + (fn*l3*len*cos(phi)^2)/(2*sin(phi)^2);
             0, 0, 0, 0];
        
        % Matrix B
        B = [l2*sin(theta), l1*cos(theta), 0;
            -l2*cos(theta), l1*sin(theta), 0;
            -(l3*len*cos(phi))/(2*sin(phi)), -(l3*len)/2, 0;
            0, 0, 1];

    elseif phi_val >= 3*pi/4 && phi_val < 5*pi/4
        % Matrix A
        A = [0, 0, - ft*l2*cos(theta) - fn*l1*sin(theta), 0;
             0, 0, fn*l1*cos(theta) - ft*l2*sin(theta), 0;
             0, 0, 0, (fn*l3*len)/2 + (fn*l3*len*sin(phi)^2)/(2*cos(phi)^2);
             0, 0, 0, 0];
        
        % Matrix B
        B = [l1*cos(theta), -l2*sin(theta), 0;
             l1*sin(theta), l2*cos(theta), 0;
             (l3*len*sin(phi))/(2*cos(phi)), -(l3*len)/2, 0;
             0, 0, 1];

    elseif phi_val >= 5*pi/4 && phi_val < 7*pi/4
        % Matrix A
        A = [0, 0, ft*l1*sin(theta) - fn*l2*cos(theta), 0;
             0, 0, - ft*l1*cos(theta) - fn*l2*sin(theta), 0;
             0, 0, 0, (fn*l3*len)/2 + (fn*l3*len*cos(phi)^2)/(2*sin(phi)^2);
             0, 0, 0, 0];
        
        % Matrix B
        B = [-l2*sin(theta), -l1*cos(theta), 0;
             l2*cos(theta), -l1*sin(theta), 0;
             -(l3*len*cos(phi))/(2*sin(phi)), -(l3*len)/2, 0;
             0, 0, 1];
    end
end