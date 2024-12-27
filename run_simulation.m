function metrics = run_simulation(Q, QN, R, mpc_timestep, duration, timestep, x_star, u_star, mu, L, radius, len, N, mass, x_initial, I_object)

    % System's parameters initialization
    num_steps = floor(duration / timestep);
    x = zeros(4, num_steps);
    x_dot = zeros(4, num_steps);
    x_ddot = zeros(3, num_steps);
    dx = zeros(4, num_steps);
    dp = zeros(3, num_steps);
    solver_times = zeros(1, num_steps);
    wrench = zeros(3, num_steps);
    time = zeros(1, num_steps);
    u = [0; 0; 0];
    x(:,1) = x_initial;
    x_start = [];

    ground_friction = zeros(3,1);

    % Define the extension for x_star and u_star
    x_star_extension = repmat(x_star(:, end), 1, N); % Repeat last column of x_star N times
    u_star_extension = zeros(size(u_star, 1), N);    % Create zero matrix for u_star
    % Append the extensions to x_star and u_star
    x_star = [x_star, x_star_extension];
    u_star = [u_star, u_star_extension];
    is_start = 1;

    for i = 1:floor(duration/timestep)
        dx(:,i) = x(:,i) - x_star(:,i);
        
        if (mod(i*timestep,mpc_timestep) == 0) || i == 1
            tic;
            try
                mpc_output = solve_MPC_MIQP(x_star(:, i:i+N), u_star(:, i:i+N), dx(:,i), mu, L, ...
                                           radius, len, N, mpc_timestep,  Q, QN, R, x_start, is_start);
            catch
                fprintf('MPC solver failed at step %d\n', i);
                mpc_output = previous_mpc_output; % Use last successful result
            end
            solver_times(i) = toc;
            is_start = 0;
            fprintf('Time is: %g\n', i * timestep);
            du(:,i) = mpc_output(4*(N+1)+1 : 4*(N+1)+3);
            z(:,i) = mpc_output(4*(N+1)+3*N+1 : 4*(N+1)+3*N+3);
            u(:,i) = du(:,i) + u_star(:,i);
            x_start = mpc_output;
        else
            u(:,i) = u(:,i-1);
            z(:,i) = z(:,i-1);
        end
    
        % Calculate parameters for the motion equation
        w = calculate_motion_model_parameters(u(:,i), x(3,i), len, radius, x(4,i));
        ground_friction_parameter = 1;
        [gr_frict, ~] = calculate_friction_with_ground(L, dp(:,i), ground_friction_parameter);
    
        ground_friction(:,i) = gr_frict;
        wrench(:,i) = w;
    
        x_dot(4,i+1) = u(3,i);
        x(4, i+1) = x(4, i) + x_dot(4, i+1) * timestep;
    
        % x_ddot(1:3, i+1) = inv(diag([mass mass I_object(3,3)])) * (-gr_frict + w);
        x_ddot(1:3, i+1) = diag([mass mass I_object(3,3)]) \ (-gr_frict + w);
        x_dot(1:3, i+1) = x_dot(1:3, i) + x_ddot(1:3, i+1) .* timestep;
        x(1:3, i+1) = x(1:3, i) + x_dot(1:3, i+1) .* timestep;
    
        time(i+1) = time(i) + timestep;
        dp(:,i+1) = [x_dot(1, i+1); x_dot(2, i+1); x_dot(3, i+1)];
    
    end

    % Calculate metrics at the end of the simulation
    metrics = calculate_metrics(x, x_star, u, solver_times, N);
end
