function [x_star_mpc, u_star_mpc] = create_mpc_star_input(x_star, u_star, ...
                                    N, iteration, timestep_parameter)
    for j = 1:N+1
        x_star_mpc(:,j) = x_star(:, iteration + (j-1)*timestep_parameter);
        u_star_mpc(:,j) = u_star(:, iteration + (j-1)*timestep_parameter);
    end
end