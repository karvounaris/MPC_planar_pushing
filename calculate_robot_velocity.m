function [v_pc, v_pc_world] = calculate_robot_velocity(L, x, len, radius, u)

    % Calculate wrench from all the contact points to the object relative to F_a
    [x_c, y_c, ~, n_c, t_c] = calculate_r_c(x(4), len, radius);
    
    J_c =  [1 0 -y_c;
            0 1 x_c];
    N = J_c' * n_c;
    T = J_c' * t_c;
    
    B = [N T];

    r_c_partial_derivative_phi = calculate_r_c_derivatives(x(4), radius);

    G_c = [J_c*L*B, r_c_partial_derivative_phi];
    v_pc = G_c * u;
    v_pc_world = [cos(x(3)) -sin(x(3)); sin(x(3)) cos(x(3))] * v_pc;

end