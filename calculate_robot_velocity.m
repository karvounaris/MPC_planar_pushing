function v_pc = calculate_robot_velocity(L, x, len, radius, u, theta_dot)

    % Calculate wrench from all the contact points to the object relative to
    % F_a
    [x_c, y_c, ~, n_c, t_c] = calculate_r_c(x(4), len, radius);
    
    J_c =  [1 0 -y_c;
            0 1 x_c];
    N = J_c' * n_c;
    T = J_c' * t_c;
    
    B = [N T];

    r_c_partial_deriative_phi = calculate_r_c_derivatives(x(4), radius);

    G_c = [J_c*L*B, r_c_partial_deriative_phi];

    v_pc = G_c * u;

    v_pc = [v_pc; 0; 0; 0; theta_dot];

    Gamma = [cos(x(3)) -sin(x(3)) 0 0 0 x(2);
             sin(x(3)) cos(x(3)) 0 0 0 -x(1);
             0 0 1 (-x(2)*cos(x(3))+x(1)*sin(x(3))) (x(2)*sin(x(3))+x(1)*cos(x(3))) 0;
             0 0 0 cos(x(3)) -sin(x(3)) 0;
             0 0 0 sin(x(3)) cos(x(3)) 0;
             0 0 0 0 0 1];

    v_pc = Gamma * v_pc;

    v_pc = [v_pc(1); v_pc(2); v_pc(6)];

end