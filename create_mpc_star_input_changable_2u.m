function [x_star_mpc, u_star_mpc, dx_mpc] = create_mpc_star_input_changable_2u(x_star, u_star, ...
                                    N, iteration, timestep_parameter, control_frequency,...
                                    u, x, len, radius, dp, timestep, L, mass, I_object)
    dx_mpc = [0; 0; 0];
    if iteration ~= 1
        for i = 1:control_frequency/timestep
            if mod(i*timestep, control_frequency) < (control_frequency/2)+1
                u_simulated = u(:,1);
            else
                u_simulated = u(:,2);
            end
            w = calculate_motion_model_parameters(u_simulated, x(3,i), len, radius, x(4,i));
            ground_friction_parameter = 1;
            [gr_frict, ~] = calculate_friction_with_ground(L, dp(:,i), ground_friction_parameter);
        
            ground_friction(:,i) = -gr_frict;
            wrench(:,i) = w;
        
            x_dot(4,i+1) = u_simulated(3);
            x(4, i+1) = x(4, i) + x_dot(4, i+1) * timestep;
        
            % x_ddot(1:3, i+1) = inv(diag([mass mass I_object(3,3)])) * (-gr_frict + w);
            x_ddot(1:3, i+1) = diag([mass mass I_object(3,3)]) \ (-gr_frict + w);
            x_dot(1:3, i+1) = x_dot(1:3, i) + x_ddot(1:3, i+1) .* timestep;
            x(1:3, i+1) = x(1:3, i) + x_dot(1:3, i+1) .* timestep;
    
            dp(:,i+1) = [x_dot(1, i+1); x_dot(2, i+1); x_dot(3, i+1)];
        end
    
        iteration = iteration + control_frequency/timestep;
        dx_mpc = x(:, end) - x_star(:, iteration);
    end

    for j = 1:N+1
        x_star_mpc(:,j) = x_star(:, iteration + (j-1)*timestep_parameter);
        u_star_mpc(:,j) = u_star(:, iteration + (j-1)*timestep_parameter);
    end
end