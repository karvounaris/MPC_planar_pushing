% % Define the MIQP solver function with Gurobi optimization
% function [mpc_output, gurobi_solve_time] = solve_MPC_MIQP(x_star, u_star, dx_0, mu, L, radius, len, N, timestep, Q, QN, R, x_start, is_start)
%     % Q = 10 * diag([3, 3, 0.1, 0]);    % State cost matrix
%     % QN = 2000 * diag([3, 3, 0.1, 0]);  % Terminal state cost matrix
%     % R = 0.5 * diag([1, 1, 0]);        % Input cost matrix
%     h = timestep;                     % Timestep
%     M = 1e2;                          % Big-M constant
% 
%     % Objective and constraints arrays
%     Aeq = [];                         % Equality constraint matrix
%     beq = [];                         % Equality constraint RHS
%     Aineq = [];                       % Inequality constraint matrix
%     bineq = [];                       % Inequality constraint RHS
% 
%     Aeq_init = [eye(4), zeros(4, 4*N + 3*N + 3*N)];
%     beq_init = dx_0;
%     Aeq = [Aeq; Aeq_init];
%     beq = [beq; beq_init];
% 
%     % Loop over horizon N to build cost and constraints
%     for i = 1:N
% 
%         % Linearize the system to get A, B at step i
%         [A, B] = linearize_system(x_star(:, i), u_star(:, i), L, radius, len);
% 
%         % Binary constraint z1 + z2 + z3 = 1
%         Aeq_binary = [zeros(1, 4*(N+1) + 3*N), zeros(1, 3*(i-1)), 1, 1, 1, zeros(1, 3*(N-i))];
%         beq_binary = 1;
%         Aeq = [Aeq; Aeq_binary];
%         beq = [beq; beq_binary];
% 
%         % Dynamics constraints at each step
%         sys_c = ss(A, B, eye(size(A)), zeros(size(B))); % Continuous system
%         sys_d = c2d(sys_c, h, 'zoh');                   % Discretize with ZOH
%         A_d = sys_d.A;
%         B_d = sys_d.B;
%         Aeq_dyn(1:4, 1:4*(N+1)+3*N+3*N) = [zeros(4, (i-1)*4), A_d -eye(4), zeros(4, (N-i)*4), ... 
%                                              zeros(4, (i-1)*3), B_d, zeros(4, (N-i)*3), ...
%                                              zeros(4, 3*N)];
%         % Dynamics constraints at each step
%         % Aeq_dyn(1:4, 1:4*(N+1)+3*N+3*N) = [zeros(4, (i-1)*4), eye(4) + h * A, -eye(4), zeros(4, (N-i)*4), ... 
%         %                                      zeros(4, (i-1)*3), h * B, zeros(4, (N-i)*3), ...
%         %                                      zeros(4, 3*N)];
%         beq_dyn = zeros(4, 1);
% 
%         Aeq = [Aeq; Aeq_dyn];
%         beq = [beq; beq_dyn];
% 
%         % Friction cone constraints
%         % Constraint 1: |f_t| <= mu * f_n
%         Aineq_friction_cone1 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -mu, 1, 0, zeros(1, (N-i)*3), ...
%                                 zeros(1, 3*N)];  % Ensure f_t <= mu * f_n
%         bineq_friction_cone1 = mu*u_star(1,i) - u_star(2,i);
%         Aineq_friction_cone2 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -mu, -1, 0, zeros(1, (N-i)*3), ...
%                                 zeros(1, 3*N)];  % Ensure f_t >= -mu * f_n
%         bineq_friction_cone2 = mu*u_star(1,i) + u_star(2,i);
% 
%         % Constraint 2: f_n >= 0
%         Aineq_friction_cone3 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -1, 0, 0, zeros(1, (N-i)*3), ...
%                                 zeros(1, 3*N)];  % Ensure f_n >= 0
%         bineq_friction_cone3 = u_star(1,i);
% 
%         % Mode constraints (Big-M constraints for mode switching)
%         % Sticking Mode (z1 = 1): phi_dot = 0
%         Aineq_stick_1 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), 0, 0, 1, zeros(1, (N-i)*3), ...
%                          zeros(1, 3*(i-1)), M, 0, 0, zeros(1, 3*(N-i))];
%         bineq_stick_1 = - u_star(3,i) + M;
%         Aineq_stick_2 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), 0, 0, -1, zeros(1, (N-i)*3), ...
%                          zeros(1, 3*(i-1)), M, 0, 0, zeros(1, 3*(N-i))];
%         bineq_stick_2 = u_star(3,i) + M;
% 
%         % Sliding Left Mode (z2 = 1): phi_dot < 0 and f_t = mu * f_n
%         Aineq_slide_left1 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -mu, 1, 0, zeros(1, (N-i)*3), ...
%                              zeros(1, 3*(i-1)), 0, M, 0, zeros(1, 3*(N-i))];
%         bineq_slide_left1 = mu*u_star(1,i) - u_star(2,i) + M;
%         Aineq_slide_left2 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), mu, -1, 0, zeros(1, (N-i)*3), ...
%                              zeros(1, 3*(i-1)), 0, M, 0, zeros(1, 3*(N-i))];
%         bineq_slide_left2 = -mu*u_star(1,i) + u_star(2,i) + M;
%         Aineq_slide_left3 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), 0, 0, 1, zeros(1, (N-i)*3), ...
%                              zeros(1, 3*(i-1)), 0, M, 0, zeros(1, 3*(N-i))];
%         bineq_slide_left3 = -u_star(3,i) - eps + M;
% 
%         % Sliding Right Mode (z3 = 1): phi_dot > 0 and f_t = -mu * f_n
%         Aineq_slide_right1 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), mu, 1, 0, zeros(1, (N-i)*3), ...
%                              zeros(1, 3*(i-1)), 0, 0, M, zeros(1, 3*(N-i))];
%         bineq_slide_right1 = -mu*u_star(1,i) - u_star(2,i) + M;
%         Aineq_slide_right2 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), -mu, -1, 0, zeros(1, (N-i)*3), ...
%                              zeros(1, 3*(i-1)), 0, 0, M, zeros(1, 3*(N-i))];
%         bineq_slide_right2 = mu*u_star(1,i) + u_star(2,i) + M;
%         Aineq_slide_right3 = [zeros(1,4*(N+1)), zeros(1, (i-1)*3), 0, 0, -1, zeros(1, (N-i)*3), ...
%                              zeros(1, 3*(i-1)), 0, 0, M, zeros(1, 3*(N-i))];
%         bineq_slide_right3 = u_star(3,i) - eps + M;
% 
%         % Add mode constraints to inequality constraints
%         Aineq = [Aineq; Aineq_friction_cone1; Aineq_friction_cone2; Aineq_friction_cone3; Aineq_stick_1; ...
%                  Aineq_stick_2; Aineq_slide_left1; Aineq_slide_left2; Aineq_slide_left3; ...
%                  Aineq_slide_right1; Aineq_slide_right2; Aineq_slide_right3];
%         bineq = [bineq; bineq_friction_cone1; bineq_friction_cone2; bineq_friction_cone3; ...
%                  bineq_stick_1; bineq_stick_2; bineq_slide_left1; bineq_slide_left2; bineq_slide_left3; ...
%                  bineq_slide_right1; bineq_slide_right2; bineq_slide_right3];
% 
%     end
% 
%     % Initialize block diagonal components for model.Q
%     Q_blocks = repmat({Q}, 1, N-1);  % Q for each time step except the terminal state
%     R_blocks = repmat({R}, 1, N);    % R for each time step for control input costs
%     W_blocks = cell(1, N);           % Initialize W_blocks as a cell array
% 
%     % Define W for each step and populate W_blocks
%     for i = 1:N
%         if i == 1 || i > 4*N/6
%             W_blocks{i} = 0.1 * diag([0, 0, 0]);
%         elseif i > 1 && i <= 2*N/6
%             W_blocks{i} = 0.1 * diag([3, 3, 3]);
%         elseif i > 2*N/6 && i <= 4*N/6
%             W_blocks{i} = 0.1 * diag([1, 1, 1]);
%         end
%     end
% 
%     % Initialize the Gurobi model as an empty structure
%     model = struct();
% 
%     % Define variable bounds and types for all variables in the model
%     % model.lb = [-inf * ones(4 * (N+1), 1); -100 * ones(3 * N, 1); zeros(3 * N, 1)];
%     % model.ub = [inf * ones(4 * (N+1), 1); 100 * ones(3 * N, 1); ones(3 * N, 1)];
%     % Initialize bounds for state variables
%     model.lb = -inf * ones(4 * (N+1), 1);
%     model.ub =  inf * ones(4 * (N+1), 1);
% 
%     % Define bounds for inputs (u = [fn, ft, phi_dot])
%     lb_u = [-100; -100; -10];
%     ub_u = [ 100;  100;  10];
% 
%     % Repeat input bounds for N steps
%     model.lb = [model.lb; repmat(lb_u, N, 1)];
%     model.ub = [model.ub; repmat(ub_u, N, 1)];
% 
%     % Add binary variable bounds
%     model.lb = [model.lb; zeros(3 * N, 1)];
%     model.ub = [model.ub; ones(3 * N, 1)];
%     model.vtype = [repmat('C', 1, 4 * (N+1) + 3 * N), repmat('B', 1, 3 * N)];
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
%     params.outputflag = 1;
%     % params.outputflag = 0;
%     params.IntFeasTol = 1e-9;
%     params.MIPGap = 0;
%     params.OptimalityTol = 1e-9;
%     params.FeasibilityTol = 1e-9;
%     params.Heuristics = 0.1;
%     params.TimeLimit = 300;
%     params.MIPFocus = 2;
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



% Define the MIQP solver function with Gurobi optimization
function [mpc_output, gurobi_solve_time] = solve_MPC_MIQP(x_star, u_star, dx_0, mu, L, radius, len, N, timestep, Q, QN, R, x_start, is_start)
    % Q = 10 * diag([3, 3, 0.1, 0]);    % State cost matrix
    % QN = 2000 * diag([3, 3, 0.1, 0]);  % Terminal state cost matrix
    % R = 0.5 * diag([1, 1, 0]);        % Input cost matrix
    h = timestep;                     % Timestep
    M = 50;                          % Big-M constant

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
        elseif i + 4 == N + 1
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
        [A, B] = jacobian_function_system_linearization(x_star(:, i), u_star(:, i), L, radius, len);

        % Dynamics constraints at each step
        sys_c = ss(A, B, eye(size(A)), zeros(size(B))); % Continuous system
        sys_d = c2d(sys_c, h, 'zoh');                   % Discretize with ZOH
        A_d = sys_d.A;
        B_d = sys_d.B;
        Aeq_dyn(1:4, 1:4*(N+1)+3*N+3*(N/5+1)) = [zeros(4, (i-1)*4), A_d -eye(4), zeros(4, (N-i)*4), ... 
                                             zeros(4, (i-1)*3), B_d, zeros(4, (N-i)*3), ...
                                             zeros(4, 3*(N/5+1))];
        % Dynamics constraints at each step
        % Aeq_dyn(1:4, 1:4*(N+1)+3*N+3*N) = [zeros(4, (i-1)*4), eye(4) + h * A, -eye(4), zeros(4, (N-i)*4), ... 
        %                                      zeros(4, (i-1)*3), h * B, zeros(4, (N-i)*3), ...
        %                                      zeros(4, 3*N)];
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

        % % Aditional robot velocity constraints
        % Gc_vel_lin = jacobian_function_velocity_linearization(x_star(:, i), u_star(:, i), L, radius, len);
        % Aineq_vel_lin1 = [zeros(1,4*(i-1)), 0, 0, 0, Gc_vel_lin(1,1), zeros(1, (N+1-i)*4), zeros(1, (i-1)*3), ...
        %                   Gc_vel_lin(1, 2:4), zeros(1, (N-i)*3), zeros(1, 3*(N/5+1))];
        % bineq_vel_lin1 = -Gc_vel_lin(1, :) * [x_star(4,i); u_star(1,i); u_star(2,i); u_star(3,i)] + 0.3;
        % Aineq_vel_lin2 = [zeros(1,4*(i-1)), 0, 0, 0, Gc_vel_lin(2,1), zeros(1, (N+1-i)*4), zeros(1, (i-1)*3), ...
        %                   Gc_vel_lin(2, 2:4), zeros(1, (N-i)*3), zeros(1, 3*(N/5+1))];
        % bineq_vel_lin2 = -Gc_vel_lin(2, :) * [x_star(4,i); u_star(1,i); u_star(2,i); u_star(3,i)] + 0.3;

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
    W_blocks = cell(1, (N/5+1));           % Initialize W_blocks as a cell array

    % Define W for each step and populate W_blocks
    for i = 1:(N/5+1)
        if i == 1 || i > N/5-1
            W_blocks{i} = 0.1 * diag([0, 0, 0]);
        elseif i == 2
            W_blocks{i} = 0.1 * diag([3, 3, 3]);
        elseif i > 2 && i <= N/5-1
            W_blocks{i} = 0.1 * diag([1, 1, 1]);
        end
    end

    % Initialize the Gurobi model as an empty structure
    model = struct();

    % Define variable bounds and types for all variables in the model
    % model.lb = [-inf * ones(4 * (N+1), 1); -100 * ones(3 * N, 1); zeros(3 * N, 1)];
    % model.ub = [inf * ones(4 * (N+1), 1); 100 * ones(3 * N, 1); ones(3 * N, 1)];
    % Initialize bounds for state variables
    model.lb = -inf * ones(4 * (N+1), 1);
    model.ub =  inf * ones(4 * (N+1), 1);

    % Define bounds for inputs (u = [fn, ft, phi_dot])
    lb_u = [-50; -30; -1];
    ub_u = [ 50;  30;  1];

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