function [solver, args] = create_MPC_MPCC_planner_controller_avoid_obstacle( ...
    Q, QN, R, w_eps0, k_eps, N, mpc_timestep, mu, L, radius, len, n_obstacles, object_shape, wid)

    import casadi.*

    x = SX.sym('x');
    y = SX.sym('y');
    theta = SX.sym('theta');
    phi = SX.sym('phi');
    states = [x; y; theta; phi];
    n_states = length(states);

    fn = SX.sym('fn');
    ft = SX.sym('ft');
    phi_dot_plus = SX.sym('phi_dot_plus');
    phi_dot_minus = SX.sym('phi_dot_minus');
    inputs = [fn; ft; phi_dot_plus; phi_dot_minus];
    n_inputs = length(inputs);

    eps = SX.sym('eps');
    n_eps = length(eps);
    
    X = SX.sym('X', n_states, (N+1));
    U = SX.sym('U', n_inputs, N);
    E = SX.sym('eps', n_eps, N);
    X_0 = SX.sym('X_0', n_states, 1);
    X_end = SX.sym('X_end', n_states, 1);  % init + reference state
    obstacles = SX.sym('Obs', 3, n_obstacles);   % [x_obs; y_obs; r_obs]

    args = struct;

    args.lbx(1:4:4*(N+1), 1) = -inf;        % lower bound for x
    args.ubx(1:4:4*(N+1), 1) = inf;         % upper bound for x
    args.lbx(2:4:4*(N+1), 1) = -inf;        % lower bound for y
    args.ubx(2:4:4*(N+1), 1) = inf;         % upper bound for y
    args.lbx(3:4:4*(N+1), 1) = -inf;        % lower bound for theta
    args.ubx(3:4:4*(N+1), 1) = inf;         % upper bound for theta
    args.lbx(4:4:4*(N+1), 1) = 5*pi/4;      % lower bound for phi
    args.ubx(4:4:4*(N+1), 1) = 7*pi/4;      % upper bound for phi

    args.lbx(4*(N+1)+1:4:4*(N+1)+4*N, 1) = 0;      % lower bound for fn
    args.ubx(4*(N+1)+1:4:4*(N+1)+4*N, 1) = 22;     % upper bound for fn
    args.lbx(4*(N+1)+2:4:4*(N+1)+4*N, 1) = -20;    % lower bound for ft
    args.ubx(4*(N+1)+2:4:4*(N+1)+4*N, 1) = 20;     % upper bound for ft
    args.lbx(4*(N+1)+3:4:4*(N+1)+4*N, 1) = 0;      % lower bound for phi_plus
    args.ubx(4*(N+1)+3:4:4*(N+1)+4*N, 1) = 2;      % upper bound for phi_plus
    args.lbx(4*(N+1)+4:4:4*(N+1)+4*N, 1) = 0;      % lower bound for phi_minus
    args.ubx(4*(N+1)+4:4:4*(N+1)+4*N, 1) = 2;      % upper bound for phi_minus

    args.lbx(4*(N+1)+4*N+1:1:4*(N+1)+4*N+N) = 0;   % lower bound for epsilon
    args.ubx(4*(N+1)+4*N+1:1:4*(N+1)+4*N+N) = inf; % upper bound for epsilon

    obj = 0;
    g = [];

    g = [g; X(:, 1) - X_0];
    args.lbg(1:4, 1) = 0;
    args.ubg(1:4, 1) = 0;

    % System dynamics function
    Rot = [cos(states(3)), -sin(states(3)), 0;
           sin(states(3)),  cos(states(3)), 0;
           0,            0,             1];
    [x_c, y_c, n_c, t_c] = calculate_r_c_casadi(states(4), len, radius, wid, object_shape);

    J_c = [1, 0, -y_c;
           0, 1,  x_c];
    B = [J_c'*n_c, J_c'*t_c];

    rhs = [Rot * L * B * [inputs(1); inputs(2)];
          inputs(3) - inputs(4)];

    f_dyn = Function('f_dyn', {states, inputs}, {rhs});
    
    constraint_offset = 5;
    
    for i = 1:N
        % 1. Dynamics (4 constraints)
        f_dyn_val = f_dyn(X(:,i), U(:,i));
        g = [g; X(:,i+1) - (X(:,i) + mpc_timestep * f_dyn_val)];
        args.lbg(constraint_offset:constraint_offset+3,1) = 0;
        args.ubg(constraint_offset:constraint_offset+3,1) = 0;
        constraint_offset = constraint_offset + 4;
    
        % 2. Lambda+ (1 constraint)
        lambda_plus = mu * U(1,i) - U(2,i);
        g = [g; lambda_plus]; 
        args.lbg(constraint_offset,1) = 0; 
        args.ubg(constraint_offset,1) = inf;      
        constraint_offset = constraint_offset + 1; 
    
        % 3. Lambda- (1 constraint)
        lambda_minus = mu * U(1,i) + U(2,i);
        g = [g; lambda_minus];
        args.lbg(constraint_offset,1) = 0;
        args.ubg(constraint_offset,1) = inf;
        constraint_offset = constraint_offset + 1;
    
        % 4. Complementarity constraint (1 constraint)
        g = [g; lambda_minus * U(3,i) + lambda_plus * U(4,i) + E(i)];
        args.lbg(constraint_offset,1) = 0;
        args.ubg(constraint_offset,1) = 0;
        constraint_offset = constraint_offset + 1;

        % obstacle constraint
        for j = 1:n_obstacles
            if j == 1
                obs_dist = sqrt((X(1,i) - obstacles(1,j))^2 + (X(2,i) - obstacles(2,j))^2);
                g = [g; obs_dist - obstacles(3,j) - len/2 - radius - 0.002];
                args.lbg(constraint_offset,1) = 0;
                args.ubg(constraint_offset,1) = inf;
                constraint_offset = constraint_offset + 1;
            end
            obs_dist = sqrt((X(1,i) - obstacles(1,j))^2 + (X(2,i) - obstacles(2,j))^2);
            g = [g; obs_dist - obstacles(3,j) - len/2 - radius - 0.002];
            args.lbg(constraint_offset,1) = 0;
            args.ubg(constraint_offset,1) = inf;
            constraint_offset = constraint_offset + 1;
        end

        % 5. Constraints for the robots velocity
        vx_object_frame = U(1,i)*(-L(3,3)*radius^2*cos(X(4,i))/sin(X(4,i))) + U(2,i)*(-L(1,1)-L(3,3)*radius^2) + (U(3,i)-U(4,i))*(radius*csc(X(4,i))^2);
        vy_object_frame = U(1,i)*(L(2,2)+L(3,3)*radius^2*(cos(X(4,i))/sin(X(4,i)))^2) + U(2,i)*(L(3,3)*radius^2*cos(X(4,i))/sin(X(4,i)));
        vx_world_frame = vx_object_frame*cos(X(3,i)) - vy_object_frame*sin(X(3,i));
        vy_world_frame = vx_object_frame*sin(X(3,i)) + vy_object_frame*cos(X(3,i));
        g = [g; vx_object_frame];
        args.lbg(constraint_offset,1) = -0.1;
        args.ubg(constraint_offset,1) = 0.1;
        constraint_offset = constraint_offset + 1;
        g = [g; vy_object_frame];
        args.lbg(constraint_offset,1) = 0;
        args.ubg(constraint_offset,1) = 0.05;
        constraint_offset = constraint_offset + 1;

        % % 5. Constraints for the object velocity
        % g = [g; (X(1,i+1) - X(1,i))/mpc_timestep; (X(2,i+1) - X(2,i))/mpc_timestep; (X(3,i+1) - X(3,i))/mpc_timestep];
        % args.lbg(constraint_offset:constraint_offset+1,1) = -0.06;
        % args.ubg(constraint_offset:constraint_offset+1,1) = 0.06;
        % constraint_offset = constraint_offset + 2;
        % args.lbg(constraint_offset,1) = -0.5;
        % args.ubg(constraint_offset,1) = 0.5;
        % constraint_offset = constraint_offset + 1;

        % 5. Constraints for the object velocity
        % vx = (X(1,i+1) - X(1,i)) / mpc_timestep;
        % vy = (X(2,i+1) - X(2,i)) / mpc_timestep;
        % v_sq = vx^2 + vy^2;
        % g = [g; v_sq];
        % args.lbg(constraint_offset,1) = 0;
        % args.ubg(constraint_offset,1) = 0.05^2;
        % constraint_offset = constraint_offset + 1;
 
    end

    for i = 1:N
        state_error = X(:,i+1) - X_end;
        % if i < N
        %     obj = obj + state_error'*Q*state_error + U(:,i)'*R*U(:,i) + w_eps0 * E(i)^2;
        % else
        %     obj = obj + state_error'*QN*state_error;
        % end
        if i < N
            obj = obj + state_error'*Q*state_error + U(:,i)'*R*U(:,i) + w_eps0 * exp(-k_eps*(i-1)) * E(i)^2;
        else
            obj = obj + state_error'*QN*state_error;
        end

        % obj = obj + state_error'*Q*state_error + U(:,i)'*R*U(:,i) + w_eps0 * exp(-k_eps*(i-1)) * E(i)^2;

        % input_stabilization = U(:,i+1) - U(:,i);
        % obj = obj + state_error'*Q*state_error + U(:,i)'*R*U(:,i) + input_stabilization'*0.1*R*input_stabilization + w_eps0 * exp(-k_eps*(i-1)) * E(i)^2;
        % obj = obj + state_error'*Q*state_error + input_stabilization'*R*input_stabilization + w_eps0 * exp(-k_eps*(i-1)) * E(i)^2;

    end

    opt_var = [reshape(X, 4*(N+1),1); reshape(U, 4*N,1); reshape(E, N,1)];

    opts = struct;
    opts.print_time = false;
    opts.verbose = false;
    opts.record_time = true;
    opts.knitro.algorithm = 2;
    opts.knitro.outlev = 0;
    opts.knitro.gradopt = 1;
    opts.knitro.hessopt = 1;
    opts.knitro.maxit = 10000;
    opts.knitro.xtol = 1e-6;
    opts.knitro.feastol = 1e-3;
    opts.knitro.feastol_abs = 1e-3;
    opts.knitro.opttol = 1e-3;
    opts.knitro.opttol_abs = 1e-3;
    opts.knitro.bar_maxcrossit = 0;
    opts.knitro.strat_warm_start = 1;
    opts.knitro.numthreads = 8;
    opts.knitro.restarts = 0;
    opts.knitro.restarts_maxit = 800;
    opts.knitro.scale = 3;
    opts.knitro.maxtime = 1;
    % opts.knitro.honorbnds = 1;

    nlp_prob = struct('f', obj, 'x', opt_var, 'g', g, 'p', [X_0(:); X_end(:); obstacles(:)]);

    solver = nlpsol('solver', 'knitro', nlp_prob, opts);
  
end

