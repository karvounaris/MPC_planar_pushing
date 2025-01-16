function v_pc = calculate_robot_velocity(L, phi, len, radius, u, theta)

    % Calculate wrench from all the contact points to the object relative to
    % F_a
    [x_c, y_c, ~, n_c, t_c] = calculate_r_c(phi, len, radius);
    
    J_c =  [1 0 -y_c;
            0 1 x_c];
    N = J_c' * n_c;
    T = J_c' * t_c;
    
    B = [N T];

    r_c_partial_deriative_phi = calculate_r_c_derivatives(phi, radius);

    G_c = [J_c*L*B, r_c_partial_deriative_phi];

    v_pc = G_c * u;

    Gamma = [cos(theta) -sin(theta);
             sin(theta) cos(theta)];

    v_pc = Gamma * v_pc;

end