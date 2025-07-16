function [solver_controller, args] = create_MPC_MIQP_controller( ...
    Q, QN, R, W, N, mpc_timestep, mu, L, radius, len, wid, object_shape, n_obstacles)

    import casadi.*

    M = 100;
    n_z_modes = 4;

    x = SX.sym('x');
    y = SX.sym('y');
    theta = SX.sym('theta');
    phi = SX.sym('phi');
    states = [x; y; theta; phi];
    n_states = length(states);

    fn = SX.sym('fn');
    ft = SX.sym('ft');
    phi_dot = SX.sym('phi_dot');
    inputs = [fn; ft; phi_dot];
    n_inputs = length(inputs);

    z1 = SX.sym('z1');
    z2 = SX.sym('z2');
    z3 = SX.sym('z3');
    z = [z1; z2; z3];
    n_z = length(z);

    % Symbolic Variables
    X = MX.sym('X', n_states, N+1);
    U = MX.sym('U', n_inputs, N);
    Z = MX.sym('Z', n_z, n_z_modes);

    X_0 = MX.sym('X_0', n_states, 1);
    X_star = MX.sym('X_star', n_states, N+1);
    U_star = MX.sym('U_star', n_inputs, N);
    obstacles = MX.sym('Obs', 3, n_obstacles);

    args = struct;

    args.lbx(1:4:4*(N+1), 1) = -inf;        % lower bound for x
    args.ubx(1:4:4*(N+1), 1) = inf;         % upper bound for x
    args.lbx(2:4:4*(N+1), 1) = -inf;        % lower bound for y
    args.ubx(2:4:4*(N+1), 1) = inf;         % upper bound for y
    args.lbx(3:4:4*(N+1), 1) = -inf;        % lower bound for theta
    args.ubx(3:4:4*(N+1), 1) = inf;         % upper bound for theta
    args.lbx(4:4:4*(N+1), 1) = 5*pi/4;      % lower bound for phi
    args.ubx(4:4:4*(N+1), 1) = 7*pi/4;      % upper bound for phi

    args.lbx(4*(N+1)+1:3:4*(N+1)+3*N, 1) = 0;      % lower bound for fn
    args.ubx(4*(N+1)+1:3:4*(N+1)+3*N, 1) = 22;     % upper bound for fn
    args.lbx(4*(N+1)+2:3:4*(N+1)+3*N, 1) = -20;    % lower bound for ft
    args.ubx(4*(N+1)+2:3:4*(N+1)+3*N, 1) = 20;     % upper bound for ft
    args.lbx(4*(N+1)+3:3:4*(N+1)+3*N, 1) = -1;      % lower bound for phi_dot
    args.ubx(4*(N+1)+3:3:4*(N+1)+3*N, 1) = 1;      % upper bound for phi_dot

    args.lbx(4*(N+1)+3*N+1 : 4*(N+1)+3*N+3*4) = 0;
    args.ubx(4*(N+1)+3*N+1 : 4*(N+1)+3*N+3*4) = 1;

    obj = 0;
    g = [];
    args.lbg = [];
    args.ubg = [];

    % Initial state constraint (equality)
    g = [g; X(:, 1) - X_0];
    args.lbg = [args.lbg; zeros(n_states, 1)];
    args.ubg = [args.ubg; zeros(n_states, 1)];

    [A_fun, B_fun] = casadi_linearize_system(L, radius, len, wid, object_shape);

    j = 0;
    for i = 1:N
        
        if i == 1 || i == 2 || i == 7 || i == 12
            j = j + 1;
            % Mode constraint: z1 + z2 + z3 = 1 (equality)
            g = [g; Z(1,j) + Z(2,j) + Z(3,j)];
            args.lbg = [args.lbg; 1];
            args.ubg = [args.ubg; 1];
        end
        
        % Dynamics constraint (equality)
        A_i = A_fun(X_star(:,i), U_star(:,i));
        B_i = B_fun(X_star(:,i), U_star(:,i));
        dyn = (X(:,i+1)-X_star(:,i+1)) - (X(:,i)-X_star(:,i)) - mpc_timestep*(A_i * (X(:,i)-X_star(:,i)) + B_i * (U(:,i)-U_star(:,i)));
        g = [g; dyn];
        args.lbg = [args.lbg; zeros(n_states,1)];
        args.ubg = [args.ubg; zeros(n_states,1)];

        % Friction cone (inequality)
        g = [g;
            -U(2,i) + mu * U(1,i);
             U(2,i) + mu * U(1,i)];
        args.lbg = [args.lbg; 0; 0];
        args.ubg = [args.ubg; inf; inf];

        % Big-M (inequality)
        g = [g;
            -U(3,i) + M*(1-Z(1,j));             % Sticking Mode (z1 = 1): phi_dot = 0
             U(3,i) + M*(1-Z(1,j));
            -U(2,i) + mu*U(1,i) + M*(1-Z(2,j)); % Sliding Left Mode (z2 = 1): phi_dot < 0 and f_t = mu * f_n
             U(2,i) - mu*U(1,i) + M*(1-Z(2,j));
            -U(3,i) + M*(1-Z(2,j)) - eps;
             U(2,i) + mu*U(1,i) + M*(1-Z(3,j)); % Sliding Right Mode (z3 = 1): phi_dot > 0 and f_t = -mu * f_n
            -U(2,i) - mu*U(1,i) + M*(1-Z(3,j));
             U(3,i) + M*(1-Z(3,j)) - eps];
        args.lbg = [args.lbg; zeros(8,1)];
        args.ubg = [args.ubg; inf(8,1)];

        % Obstacle avoidance
        for k = 1:n_obstacles
            obs_dist = sqrt((X(1,i) - obstacles(1,k))^2 + (X(2,i) - obstacles(2,k))^2);
            g = [g; obs_dist - obstacles(3,k) - len/2 - radius - 0.001];
            args.lbg = [args.lbg; 0];
            args.ubg = [args.ubg; inf];
        end
        
        % Cost
        state_error = X(:,i+1) - X_star(:, i+1);
        input_error = U(:,i) - U_star(:,i);
        if i < N
            obj = obj + state_error'*Q*state_error + input_error'*R*input_error + Z(:,j)'*W*Z(:,j);
        else
            obj = obj + state_error'*QN*state_error;
        end
    end
        
    opt_var = [reshape(X, n_states*(N+1), 1); ...
               reshape(U, n_inputs*N, 1); ...
               reshape(Z, n_z*n_z_modes, 1)];

    param = [reshape(X_0, n_states, 1); ...
             reshape(X_star, n_states*(N+1), 1); ...
             reshape(U_star, n_inputs*N, 1); ...
             reshape(obstacles, 3*n_obstacles, 1)];

    bin_var_indices = 4*(N+1) + 3*N + (1:3*n_z_modes);
    n_opt = length(opt_var);

    qp = struct;
    qp.x = opt_var;
    qp.p = param;
    qp.f = obj;
    qp.g = g; 

    qp_opts = struct;
    qp_opts.print_time = false;
    qp_opts.verbose = true;
    qp_opts.record_time = true;
    qp_opts.gurobi.TimeLimit = 10;
    qp_opts.gurobi.OutputFlag = 1;
    qp_opts.gurobi.Threads = 6;
    qp_opts.discrete = false(n_opt, 1);
    qp_opts.discrete(bin_var_indices) = true;

    solver_controller = qpsol('solver', 'gurobi', qp, qp_opts);

end

