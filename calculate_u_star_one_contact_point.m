%=======================================================================%
% This function solves a 3x3 system and finds the u_star, given:
%     - L, matrix of limit surface
%     - x_star
%     - y_star
%     - theta_star
%     - len, length of object (for more check calculate_r_c.m)
%     - radius, radius of object (for more check calculate_r_c.m)
%=======================================================================%

function [fn_star, ft_star, phi_star_dot, phi_star] = calculate_u_star_one_contact_point(L, x_star_dot, y_star_dot, theta_star, ...
                                                        theta_star_dot, len, radius, timestep, duration, wid, object_shape)
    

    options = optimoptions('fsolve', ...
        'MaxFunctionEvaluations', 1e6, ...
        'MaxIterations', 4e6, ...
        'FiniteDifferenceType', 'central', ...   % 'central' or 'forward'
        'FunctionTolerance', 1e-6, ...
        'StepTolerance', 1e-6, ...
        'FiniteDifferenceStepSize', 1e-4, ...
        'Algorithm', 'trust-region-dogleg', ...  % 'trust-region-dogleg' or 'levenberg-marquardt'
        'Diagnostics', 'off');
    
    for i = 1:(duration/timestep)+1
        initial_guess = [5, 0.1, 1.5];
        [solution, fval, exitflag, output] = fsolve(@(u_star) system_equations ... 
                                            (u_star, L, x_star_dot(i), y_star_dot(i), ...
                                            theta_star(i), theta_star_dot(i), len, radius, wid, object_shape), ...
                                            initial_guess, options);
        disp(output);

        fn_star(i) = exp(solution(1));
        ft_star(i) = solution(2);
        phi_star(i) = exp(solution(3));
        phi_star(i) = mod(phi_star(i), 2*pi);
    end

    for i = 1:length(phi_star)-1
        phi_star_dot(i) = (phi_star(i+1) - phi_star(i))/timestep;
    end
    phi_star_dot(i+1) = 0;
end

%=======================================================================%
% This function is the formula of the 3x3 of the system equation, the 
% of which provides the u_star, given:
%     - L, matrix of limit surface
%     - x_star
%     - y_star
%     - theta_star
%     - len, length of object (for more check calculate_r_c.m)
%     - radius, radius of object (for more check calculate_r_c.m)
%=======================================================================%

function output = system_equations(u_star, L, x_star_dot, y_star_dot, theta_star, theta_star_dot, len, radius, wid, object_shape)

    z1 = u_star(1);
    f_t_star = u_star(2);
    z2 = u_star(3);

    f_n_star = exp(z1);
    phi_star = exp(z2);

    [x_c, y_c, r_c, n_c, t_c] = calculate_r_c(phi_star, len, radius, wid, object_shape);

    a = n_c(1);
    b = n_c(2);
    c = t_c(1);
    d = t_c(2);

    output(1) = (L(1,1)*cos(theta_star)*a - L(2,2)*sin(theta_star)*b)*f_n_star ...
                + (L(1,1)*cos(theta_star)*c - L(2,2)*sin(theta_star)*d)*f_t_star - x_star_dot;
    output(2) = (L(1,1)*sin(theta_star)*a + L(2,2)*cos(theta_star)*b)*f_n_star ...
                + (L(1,1)*sin(theta_star)*c + L(2,2)*cos(theta_star)*d)*f_t_star - y_star_dot;
    output(3) = (L(3,3)*(-y_c*a+x_c*b))*f_n_star + (L(3,3)*(-y_c*c+x_c*d))*f_t_star - theta_star_dot;
end