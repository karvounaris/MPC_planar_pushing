% % Define the MIQP solver function with Gurobi optimization
% function [mpc_output, gurobi_solve_time] = solve_MPC_MIQP(x_star, u_star, dx_0, mu, L, radius, ...
%                                                           len, wid, N, timestep, Q, QN, R, W, ...
%                                                           x_start, is_start, object_shape)
%     h = timestep;                     % Timestep
%     M = 100;                          % Big-M constant
% 
%     % Objective and constraints arrays
%     Aeq = [];                         % Equality constraint matrix
%     beq = [];                         % Equality constraint RHS
%     Aineq = [];                       % Inequality constraint matrix
%     bineq = [];                       % Inequality constraint RHS
% 
%     Aeq_init = [eye(4), zeros(4, 4*N + 3*N + 3*4)];
%     beq_init = dx_0;
%     Aeq = [Aeq; Aeq_init];
%     beq = [beq; beq_init];
% 
%     j = 0;
%     % Loop over horizon N to build cost and constraints
%     for i = 1:N
%         tic
%         if i == 1 || i == 2 || i == 7 || i == 12
%             j = j + 1;
%             % Binary constraint z1 + z2 + z3 = 1
%             Aeq_binary = [zeros(1, 4*(N+1) + 3*N), zeros(1, 3*(j-1)), 1, 1, 1, zeros(1, 3*(4-j))];
%             beq_binary = 1;
%             Aeq = [Aeq; Aeq_binary];
%             beq = [beq; beq_binary];
% 
%         % [Ai, Bi] = linearize_system(x_star(:, i), u_star(:, i), L, radius, len);
%         [A, B] = jacobian_function_system_linearization(x_star(:, i), u_star(:, i), L, radius, len, wid, object_shape);
% 
%         % Dynamics constraints at each step
%         Aeq_dyn(1:4, 1:4*(N+1)+3*N+3*4) = [zeros(4, (i-1)*4), eye(4) + h * A, -eye(4), zeros(4, (N-i)*4), ... 
%                                              zeros(4, (i-1)*3), h * B, zeros(4, (N-i)*3), ...
%                                              zeros(4, 3*4)];
%         beq_dyn = zeros(4, 1);
% 
%         Aeq = [Aeq; Aeq_dyn];
%         beq = [beq; beq_dyn];
% 
%         % Friction cone constraints
%         % Constraint 1: |f_t| <= mu * f_n
%         Aineq_friction_cone1 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -mu, 1, 0, zeros(1, (N-i)*3), ...
%                                 zeros(1, 3*4)];  % Ensure f_t <= mu * f_n
%         bineq_friction_cone1 = mu*u_star(1,i) - u_star(2,i);
%         Aineq_friction_cone2 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -mu, -1, 0, zeros(1, (N-i)*3), ...
%                                 zeros(1, 3*4)];  % Ensure f_t >= -mu * f_n
%         bineq_friction_cone2 = mu*u_star(1,i) + u_star(2,i);
% 
%         % Constraint 2: f_n >= 0
%         Aineq_friction_cone3 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -1, 0, 0, zeros(1, (N-i)*3), ...
%                                 zeros(1, 3*4)];  % Ensure f_n >= 0
%         bineq_friction_cone3 = u_star(1,i);
% 
%         % Mode constraints (Big-M constraints for mode switching)
%         % Sticking Mode (z1 = 1): phi_dot = 0
%         Aineq_stick_1 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), 0, 0, 1, zeros(1, (N-i)*3), ...
%                          zeros(1, 3*(j-1)), M, 0, 0, zeros(1, 3*(4-j))];
%         bineq_stick_1 = - u_star(3,i) + M;
%         Aineq_stick_2 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), 0, 0, -1, zeros(1, (N-i)*3), ...
%                          zeros(1, 3*(j-1)), M, 0, 0, zeros(1, 3*(4-j))];
%         bineq_stick_2 = u_star(3,i) + M;
% 
%         % Sliding Left Mode (z2 = 1): phi_dot < 0 and f_t = mu * f_n
%         Aineq_slide_left1 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -mu, 1, 0, zeros(1, (N-i)*3), ...
%                              zeros(1, 3*(j-1)), 0, M, 0, zeros(1, 3*(4-j))];
%         bineq_slide_left1 = mu*u_star(1,i) - u_star(2,i) + M;
%         Aineq_slide_left2 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), mu, -1, 0, zeros(1, (N-i)*3), ...
%                              zeros(1, 3*(j-1)), 0, M, 0, zeros(1, 3*(4-j))];
%         bineq_slide_left2 = -mu*u_star(1,i) + u_star(2,i) + M;
%         Aineq_slide_left3 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), 0, 0, 1, zeros(1, (N-i)*3), ...
%                              zeros(1, 3*(j-1)), 0, M, 0, zeros(1, 3*(4-j))];
%         bineq_slide_left3 = -u_star(3,i) - eps + M;
% 
%         % Sliding Right Mode (z3 = 1): phi_dot > 0 and f_t = -mu * f_n
%         Aineq_slide_right1 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), mu, 1, 0, zeros(1, (N-i)*3), ...
%                              zeros(1, 3*(j-1)), 0, 0, M, zeros(1, 3*(4-j))];
%         bineq_slide_right1 = -mu*u_star(1,i) - u_star(2,i) + M;
%         Aineq_slide_right2 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -mu, -1, 0, zeros(1, (N-i)*3), ...
%                              zeros(1, 3*(j-1)), 0, 0, M, zeros(1, 3*(4-j))];
%         bineq_slide_right2 = mu*u_star(1,i) + u_star(2,i) + M;
%         Aineq_slide_right3 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), 0, 0, -1, zeros(1, (N-i)*3), ...
%                              zeros(1, 3*(j-1)), 0, 0, M, zeros(1, 3*(4-j))];
%         bineq_slide_right3 = u_star(3,i) - eps + M;
% 
%         % Add mode constraints to inequality constraints
%         Aineq = [Aineq; Aineq_friction_cone1; Aineq_friction_cone2; Aineq_friction_cone3; Aineq_stick_1; ...
%                  Aineq_stick_2; Aineq_slide_left1; Aineq_slide_left2; Aineq_slide_left3; ...
%                  Aineq_slide_right1; Aineq_slide_right2; Aineq_slide_right3];
%         bineq = [bineq; bineq_friction_cone1; bineq_friction_cone2; bineq_friction_cone3; ...
%                  bineq_stick_1; bineq_stick_2; bineq_slide_left1; bineq_slide_left2; bineq_slide_left3; ...
%                  bineq_slide_right1; bineq_slide_right2; bineq_slide_right3];
%     end
% 
%     % Initialize block diagonal components for model.Q
%     Q_blocks = repmat({Q}, 1, N-1);  % Q for each time step except the terminal state
%     R_blocks = repmat({R}, 1, N);    % R for each time step for control input costs
%     W_blocks = repmat({W}, 1, 4); % Initialize W_blocks as a cell array
% 
%     % Initialize the Gurobi model as an empty structure
%     model = struct();
% 
%     % Initialize bounds for state variables
%     model.lb = -inf * ones(4 * (N+1), 1);
%     model.ub =  inf * ones(4 * (N+1), 1);
% 
%     % Define bounds for inputs (u = [fn, ft, phi_dot])
%     lb_u = [-30; -30; -3];
%     ub_u = [4;  30;  3];
% 
%     % Repeat input bounds for N steps
%     model.lb = [model.lb; repmat(lb_u, N, 1)];
%     model.ub = [model.ub; repmat(ub_u, N, 1)];
% 
%     % Add binary variable bounds
%     model.lb = [model.lb; zeros(3 * 4, 1)];
%     model.ub = [model.ub; ones(3 * 4, 1)];
%     model.vtype = [repmat('C', 1, 4 * (N+1) + 3 * N), repmat('B', 1, 3 * 4)];
% 
%     % Set the remaining fields in model (A, rhs, Q, etc.) as in your code
%     model.A = sparse([Aeq; Aineq]);
%     model.rhs = [beq; bineq];
%     model.sense = [repmat('=', size(beq, 1), 1); repmat('<', size(bineq, 1), 1)];
%     model.Q = sparse(blkdiag(zeros(4,4), Q_blocks{:}, QN, R_blocks{:}, W_blocks{:}));
%     if is_start == 0
%         model.start = x_start;
%     end
% 
%     % Additional Gurobi parameters
%     % params.outputflag = 1;
%     params.outputflag = 0;
%     % params.IntFeasTol = 1e-9;
%     % params.MIPGap = 0;
%     % params.OptimalityTol = 1e-9;
%     % params.FeasibilityTol = 1e-9;
%     % params.Heuristics = 0.1;
%     params.TimeLimit = 300;
%     % params.MIPFocus = 2;
%     params.Threads = 2;
%     % model.Params.MemLimit = 1;
% 
%     % Solve the problem with Gurobi
%     result = gurobi(model, params);
% 
%     % Retrieve the solve time from Gurobi's output
%     gurobi_solve_time = result.runtime;
% 
%     % Check the solver status and process the result
%     disp(result.status);
%     if strcmp(result.status, 'OPTIMAL')
%         % Extract the optimal solution
%         mpc_output = result.x;
%     elseif strcmp(result.status, 'TIME_LIMIT')
%         warning('Gurobi did not find an optimal solution within the time limit. Status: %s', result.status);
%         mpc_output = [];
%     else
%         error('Gurobi solver failed with status: %s', result.status);
%     end
% 
%     % Check for solution feasibility
%     if strcmp(result.status, 'OPTIMAL')
%         % Extract optimal control input
%         mpc_output = result.x;
%     elseif strcmp(result.status, 'TIME_LIMIT')
%         warning('Gurobi did not find an optimal solution within the time limit. Status: %s', result.status);
%         u_opt = [];
%         z = [];
%     else
%         error('Gurobi solver failed with status: %s', result.status);
%     end
% end
% 





































% Define the MIQP solver function with Gurobi optimization
function [mpc_output, gurobi_solve_time] = solve_MPC_MIQP(x_star, u_star, dx_0, mu, L, radius, ...
                                                          len, wid, N, timestep, Q, QN, R, W, ...
                                                          x_start, is_start, object_shape)
    h = timestep;                     % Timestep
    M = 100;                          % Big-M constant

    % Objective and constraints arrays
    Aeq = [];                         % Equality constraint matrix
    beq = [];                         % Equality constraint RHS
    Aineq = [];                       % Inequality constraint matrix
    bineq = [];                       % Inequality constraint RHS

    Aeq_init = [eye(4), zeros(4, 4*N + 3*N + 3*(N/5+1))];
    beq_init = dx_0;
    Aeq = [Aeq; Aeq_init];
    beq = [beq; beq_init];

    j = 0;
    % Loop over horizon N to build cost and constraints
    for i = 1:N
        tic
        if i == 1
            j = j + 1;
            % Binary constraint z1 + z2 + z3 = 1
            Aeq_binary = [zeros(1, 4*(N+1) + 3*N), zeros(1, 3*(j-1)), 1, 1, 1, zeros(1, 3*((N/5+1)-j))];
            beq_binary = 1;
            Aeq = [Aeq; Aeq_binary];
            beq = [beq; beq_binary];
        elseif mod(i, 5) == 2
            j = j + 1;
            % Binary constraint z1 + z2 + z3 = 1
            Aeq_binary = [zeros(1, 4*(N+1) + 3*N), zeros(1, 3*(j-1)), 1, 1, 1, zeros(1, 3*((N/5+1)-j))];
            beq_binary = 1;
            Aeq = [Aeq; Aeq_binary];
            beq = [beq; beq_binary];
        end

        % [Ai, Bi] = linearize_system(x_star(:, i), u_star(:, i), L, radius, len);
        [A, B] = jacobian_function_system_linearization(x_star(:, i), u_star(:, i), L, radius, len, wid, object_shape);

        % Dynamics constraints at each step
        Aeq_dyn(1:4, 1:4*(N+1)+3*N+3*(N/5+1)) = [zeros(4, (i-1)*4), eye(4) + h * A, -eye(4), zeros(4, (N-i)*4), ... 
                                             zeros(4, (i-1)*3), h * B, zeros(4, (N-i)*3), ...
                                             zeros(4, 3*(N/5+1))];
        beq_dyn = zeros(4, 1);

        Aeq = [Aeq; Aeq_dyn];
        beq = [beq; beq_dyn];

        % Friction cone constraints
        % Constraint 1: |f_t| <= mu * f_n
        Aineq_friction_cone1 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -mu, 1, 0, zeros(1, (N-i)*3), ...
                                zeros(1, 3*(N/5+1))];  % Ensure f_t <= mu * f_n
        bineq_friction_cone1 = mu*u_star(1,i) - u_star(2,i);
        Aineq_friction_cone2 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -mu, -1, 0, zeros(1, (N-i)*3), ...
                                zeros(1, 3*(N/5+1))];  % Ensure f_t >= -mu * f_n
        bineq_friction_cone2 = mu*u_star(1,i) + u_star(2,i);

        % Constraint 2: f_n >= 0
        Aineq_friction_cone3 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -1, 0, 0, zeros(1, (N-i)*3), ...
                                zeros(1, 3*(N/5+1))];  % Ensure f_n >= 0
        bineq_friction_cone3 = u_star(1,i);

        % Mode constraints (Big-M constraints for mode switching)
        % Sticking Mode (z1 = 1): phi_dot = 0
        Aineq_stick_1 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), 0, 0, 1, zeros(1, (N-i)*3), ...
                         zeros(1, 3*(j-1)), M, 0, 0, zeros(1, 3*((N/5+1)-j))];
        bineq_stick_1 = - u_star(3,i) + M;
        Aineq_stick_2 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), 0, 0, -1, zeros(1, (N-i)*3), ...
                         zeros(1, 3*(j-1)), M, 0, 0, zeros(1, 3*((N/5+1)-j))];
        bineq_stick_2 = u_star(3,i) + M;

        % Sliding Left Mode (z2 = 1): phi_dot < 0 and f_t = mu * f_n
        Aineq_slide_left1 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -mu, 1, 0, zeros(1, (N-i)*3), ...
                             zeros(1, 3*(j-1)), 0, M, 0, zeros(1, 3*((N/5+1)-j))];
        bineq_slide_left1 = mu*u_star(1,i) - u_star(2,i) + M;
        Aineq_slide_left2 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), mu, -1, 0, zeros(1, (N-i)*3), ...
                             zeros(1, 3*(j-1)), 0, M, 0, zeros(1, 3*((N/5+1)-j))];
        bineq_slide_left2 = -mu*u_star(1,i) + u_star(2,i) + M;
        Aineq_slide_left3 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), 0, 0, 1, zeros(1, (N-i)*3), ...
                             zeros(1, 3*(j-1)), 0, M, 0, zeros(1, 3*((N/5+1)-j))];
        bineq_slide_left3 = -u_star(3,i) - eps + M;

        % Sliding Right Mode (z3 = 1): phi_dot > 0 and f_t = -mu * f_n
        Aineq_slide_right1 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), mu, 1, 0, zeros(1, (N-i)*3), ...
                             zeros(1, 3*(j-1)), 0, 0, M, zeros(1, 3*((N/5+1)-j))];
        bineq_slide_right1 = -mu*u_star(1,i) - u_star(2,i) + M;
        Aineq_slide_right2 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -mu, -1, 0, zeros(1, (N-i)*3), ...
                             zeros(1, 3*(j-1)), 0, 0, M, zeros(1, 3*((N/5+1)-j))];
        bineq_slide_right2 = mu*u_star(1,i) + u_star(2,i) + M;
        Aineq_slide_right3 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), 0, 0, -1, zeros(1, (N-i)*3), ...
                             zeros(1, 3*(j-1)), 0, 0, M, zeros(1, 3*((N/5+1)-j))];
        bineq_slide_right3 = u_star(3,i) - eps + M;

        % Add mode constraints to inequality constraints
        Aineq = [Aineq; Aineq_friction_cone1; Aineq_friction_cone2; Aineq_friction_cone3; Aineq_stick_1; ...
                 Aineq_stick_2; Aineq_slide_left1; Aineq_slide_left2; Aineq_slide_left3; ...
                 Aineq_slide_right1; Aineq_slide_right2; Aineq_slide_right3]; ...
                 % Aineq_vel_lin1; Aineq_vel_lin2];
        bineq = [bineq; bineq_friction_cone1; bineq_friction_cone2; bineq_friction_cone3; ...
                 bineq_stick_1; bineq_stick_2; bineq_slide_left1; bineq_slide_left2; bineq_slide_left3; ...
                 bineq_slide_right1; bineq_slide_right2; bineq_slide_right3]; ...
                 % bineq_vel_lin1; bineq_vel_lin2];
    end

    % Initialize block diagonal components for model.Q
    Q_blocks = repmat({Q}, 1, N-1);  % Q for each time step except the terminal state
    R_blocks = repmat({R}, 1, N);    % R for each time step for control input costs
    W_blocks = repmat({W}, 1, (N/5+1)); % Initialize W_blocks as a cell array


    % Initialize the Gurobi model as an empty structure
    model = struct();

    % Initialize bounds for state variables
    model.lb = -inf * ones(4 * (N+1), 1);
    model.ub =  inf * ones(4 * (N+1), 1);

    % Define bounds for inputs (u = [fn, ft, phi_dot])
    lb_u = [-30; -30; -3];
    ub_u = [1;  30;  3];

    % Repeat input bounds for N steps
    model.lb = [model.lb; repmat(lb_u, N, 1)];
    model.ub = [model.ub; repmat(ub_u, N, 1)];

    % Add binary variable bounds
    model.lb = [model.lb; zeros(3 * (N/5+1), 1)];
    model.ub = [model.ub; ones(3 * (N/5+1), 1)];
    model.vtype = [repmat('C', 1, 4 * (N+1) + 3 * N), repmat('B', 1, 3 * (N/5+1))];

    % Set the remaining fields in model (A, rhs, Q, etc.) as in your code
    model.A = sparse([Aeq; Aineq]);
    model.rhs = [beq; bineq];
    model.sense = [repmat('=', size(beq, 1), 1); repmat('<', size(bineq, 1), 1)];
    model.Q = sparse(blkdiag(zeros(4,4), Q_blocks{:}, QN, R_blocks{:}, W_blocks{:}));
    if is_start == 0
        model.start = x_start;
    end

    % Additional Gurobi parameters
    % params.outputflag = 1;
    params.outputflag = 0;
    % params.IntFeasTol = 1e-9;
    % params.MIPGap = 0;
    % params.OptimalityTol = 1e-9;
    % params.FeasibilityTol = 1e-9;
    % params.Heuristics = 0.1;
    params.TimeLimit = 300;
    % params.MIPFocus = 2;
    params.Threads = 2;
    % model.Params.MemLimit = 1;

    % Solve the problem with Gurobi
    result = gurobi(model, params);

    % Retrieve the solve time from Gurobi's output
    gurobi_solve_time = result.runtime;

    % Check the solver status and process the result
    disp(result.status);
    if strcmp(result.status, 'OPTIMAL')
        % Extract the optimal solution
        mpc_output = result.x;
    elseif strcmp(result.status, 'TIME_LIMIT')
        warning('Gurobi did not find an optimal solution within the time limit. Status: %s', result.status);
        mpc_output = [];
    else
        error('Gurobi solver failed with status: %s', result.status);
    end

    % Check for solution feasibility
    if strcmp(result.status, 'OPTIMAL')
        % Extract optimal control input
        mpc_output = result.x;
    elseif strcmp(result.status, 'TIME_LIMIT')
        warning('Gurobi did not find an optimal solution within the time limit. Status: %s', result.status);
        u_opt = [];
        z = [];
    else
        error('Gurobi solver failed with status: %s', result.status);
    end
end
