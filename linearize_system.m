function [A, B] = linearize_system(x_star, u_star, L, radius, len)
    % Define symbolic variables for state and input
    syms x_sym y_sym theta_sym phi_sym fn_sym ft_sym phi_dot_sym
    x_vec = [x_sym; y_sym; theta_sym; phi_sym];
    u_vec = [fn_sym; ft_sym; phi_dot_sym];
    
    % Define symbolic expressions for each case in the dynamics
    dx_case1 = [
        (L(1,1)*(-cos(2*phi_sym))*cos(theta_sym) - L(2,2)*(-sin(2*phi_sym))*sin(theta_sym))*fn_sym + ...
        (L(1,1)*(sin(2*phi_sym))*cos(theta_sym) - L(2,2)*(-cos(2*phi_sym))*sin(theta_sym))*ft_sym;
        
        (L(1,1)*(-cos(2*phi_sym))*sin(theta_sym) + L(2,2)*(-sin(2*phi_sym))*cos(theta_sym))*fn_sym + ...
        (L(1,1)*(sin(2*phi_sym))*sin(theta_sym) + L(2,2)*(-cos(2*phi_sym))*cos(theta_sym))*ft_sym;
        
        L(3,3)*(-radius*sin(2*phi_sym)*(-cos(2*phi_sym)) + (len/2+radius*cos(2*phi_sym))*(-sin(2*phi_sym)))*fn_sym + ...
        L(3,3)*(-radius*sin(2*phi_sym)*(sin(2*phi_sym)) + (len/2+radius*cos(2*phi_sym))*(-cos(2*phi_sym)))*ft_sym;
        
        phi_dot_sym];

    dx_case2 = [
        (L(1,1)*(0)*cos(theta_sym) - L(2,2)*(-1)*sin(theta_sym))*fn_sym + ...
        (L(1,1)*(1)*cos(theta_sym) - L(2,2)*(0)*sin(theta_sym))*ft_sym;

        (L(1,1)*(0)*sin(theta_sym) + L(2,2)*(-1)*cos(theta_sym))*fn_sym + ...
        (L(1,1)*(1)*sin(theta_sym) + L(2,2)*(0)*cos(theta_sym))*ft_sym;

        L(3,3)*(-radius*(0) + radius*(cos(phi_sym)/sin(phi_sym))*(-1))*fn_sym + ...
        L(3,3)*(-radius*(1) + radius*(cos(phi_sym)/sin(phi_sym))*(0))*ft_sym;

        phi_dot_sym];

    dx_case3 = [
        (L(1,1)*(cos(2*phi_sym))*cos(theta_sym) - L(2,2)*(sin(2*phi_sym))*sin(theta_sym))*fn_sym + ...
        (L(1,1)*(-sin(2*phi_sym))*cos(theta_sym) - L(2,2)*(cos(2*phi_sym))*sin(theta_sym))*ft_sym;

        (L(1,1)*(cos(2*phi_sym))*sin(theta_sym) + L(2,2)*(sin(2*phi_sym))*cos(theta_sym))*fn_sym + ...
        (L(1,1)*(-sin(2*phi_sym))*sin(theta_sym) + L(2,2)*(cos(2*phi_sym))*cos(theta_sym))*ft_sym;

        L(3,3)*(radius*sin(2*phi_sym)*(cos(2*phi_sym)) - (len/2+radius*cos(2*phi_sym))*(sin(2*phi_sym)))*fn_sym + ...
        L(3,3)*(radius*sin(2*phi_sym)*(-sin(2*phi_sym)) - (len/2+radius*cos(2*phi_sym))*(cos(2*phi_sym)))*ft_sym;

        phi_dot_sym];

    dx_case4 = [
        (L(1,1)*(0)*cos(theta_sym) - L(2,2)*(1)*sin(theta_sym))*fn_sym + ...
        (L(1,1)*(-1)*cos(theta_sym) - L(2,2)*(0)*sin(theta_sym))*ft_sym;

        (L(1,1)*(0)*sin(theta_sym) + L(2,2)*(1)*cos(theta_sym))*fn_sym + ...
        (L(1,1)*(-1)*sin(theta_sym) + L(2,2)*(0)*cos(theta_sym))*ft_sym;

        L(3,3)*(radius*(0) - radius*(cos(phi_sym)/sin(phi_sym))*(1))*fn_sym + ...
        L(3,3)*(radius*(-1) - radius*(cos(phi_sym)/sin(phi_sym))*(0))*ft_sym;

        phi_dot_sym];

    % Select the appropriate case based on phi value
    phi_val = x_star(4);
    if (phi_val >= 0 && phi_val < pi/4) || (phi_val >= 7*pi/4 && phi_val < 2*pi)
        f_sym = dx_case1;
    elseif phi_val >= pi/4 && phi_val < 3*pi/4
        f_sym = dx_case2;
    elseif phi_val >= 3*pi/4 && phi_val < 5*pi/4
        f_sym = dx_case3;
    elseif phi_val >= 5*pi/4 && phi_val < 7*pi/4
        f_sym = dx_case4;
    end

    % Calculate the Jacobians of f_sym with respect to x and u
    A_sym = jacobian(f_sym, x_vec);
    B_sym = jacobian(f_sym, u_vec);

    % Substitute x_star and u_star into the Jacobian matrices
    A = double(subs(A_sym, {x_sym, y_sym, theta_sym, phi_sym, fn_sym, ft_sym, phi_dot_sym}, ...
                    {x_star(1), x_star(2), x_star(3), x_star(4), u_star(1), u_star(2), u_star(3)}));
    B = double(subs(B_sym, {x_sym, y_sym, theta_sym, phi_sym, fn_sym, ft_sym, phi_dot_sym}, ...
                    {x_star(1), x_star(2), x_star(3), x_star(4), u_star(1), u_star(2), u_star(3)}));
end
