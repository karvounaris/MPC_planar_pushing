function [path_error, x_y_error, theta_error] = path_error_MIQP(x, x_star)

    alpha = 0.05;
    path_error = [];
    x_y_error = [];
    theta_error = [];
    time_horizon = 250;

    for i = 1:length(x(1,:))
        
        if i <= time_horizon
            d = (x(1,i) - x_star(1,1:i+time_horizon)).^2 + (x(2,i) - x_star(2,1:i+time_horizon)).^2 + alpha * (x(3,i) - x_star(3,1:i+time_horizon)).^2;
            N = time_horizon + i;
        elseif i <= length(x(1,:)) - time_horizon
            d = (x(1,i) - x_star(1,i-time_horizon:i+time_horizon)).^2 + (x(2,i) - x_star(2,i-time_horizon:i+time_horizon)).^2 + alpha * (x(3,i) - x_star(3,i-time_horizon:i+time_horizon)).^2;
            N = 2*time_horizon + 1;
        elseif i > length(x(1,:)) - time_horizon
            d = (x(1,i) - x_star(1,i-time_horizon:length(x(1,:)))).^2 + (x(2,i) - x_star(2,i-time_horizon:length(x(1,:)))).^2 + alpha * (x(3,i) - x_star(3,i-time_horizon:length(x(1,:)))).^2;
            N = time_horizon + (length(x(1,:)) - i) + 1;
        end

        % Define Q as a diagonal matrix
        Q = diag(d);
        
        % Define constraint matrix A and vector b
        A = ones(1, N);
        b = 1;
        
        % Define binary constraints
        lb = zeros(N, 1);
        ub = ones(N, 1);
        
        % Solve MIQP with Gurobi
        model.Q = sparse(Q);
        model.A = sparse(A);
        model.rhs = b;
        model.sense = '=';
        model.lb = lb;
        model.ub = ub;
        model.vtype = 'B'; % Binary variables
        model.modelsense = 'min';
        
        % params.outputflag = 1;
        params.outputflag = 0;
        params.TimeLimit = 300;
        params.Threads = 2;

        result = gurobi(model, params);
        path_error = [path_error, sqrt(result.objval)]; % Path error

        fprintf('Set 1, Iteration %d from %d and solver time is %2.4f\n', i, length(x(1,:)), result.runtime);
    end

    for i = 1:length(x(1,:))

        if i <= time_horizon
            d = (x(1,i) - x_star(1,1:i+time_horizon)).^2 + (x(2,i) - x_star(2,1:i+time_horizon)).^2;
            N = time_horizon + i;
        elseif i <= length(x(1,:)) - time_horizon
            d = (x(1,i) - x_star(1,i-time_horizon:i+time_horizon)).^2 + (x(2,i) - x_star(2,i-time_horizon:i+time_horizon)).^2;
            N = 2*time_horizon + 1;
        elseif i > length(x(1,:)) - time_horizon
            d = (x(1,i) - x_star(1,i-time_horizon:length(x(1,:)))).^2 + (x(2,i) - x_star(2,i-time_horizon:length(x(1,:)))).^2;
            N = time_horizon + (length(x(1,:)) - i) + 1;
        end

        % Define Q as a diagonal matrix
        Q = diag(d);
        
        % Define constraint matrix A and vector b
        A = ones(1, N);
        b = 1;
        
        % Define binary constraints
        lb = zeros(N, 1);
        ub = ones(N, 1);
        
        % Solve MIQP with Gurobi
        model.Q = sparse(Q);
        model.A = sparse(A);
        model.rhs = b;
        model.sense = '=';
        model.lb = lb;
        model.ub = ub;
        model.vtype = 'B'; % Binary variables
        model.modelsense = 'min';
        
        % params.outputflag = 1;
        params.outputflag = 0;
        params.TimeLimit = 300;
        params.Threads = 2;

        result = gurobi(model, params);
        x_y_error = [x_y_error, sqrt(result.objval)]; % Path error
        fprintf('Set 2, Iteration %d from %d and solver time is %2.4f\n', i, length(x(1,:)), result.runtime);
    end

    for i = 1:length(x(1,:))
        
        if i <= time_horizon
            d = (x(3,i) - x_star(3,1:i+time_horizon)).^2;
            N = time_horizon + i;
        elseif i <= length(x(1,:)) - time_horizon
            d = (x(3,i) - x_star(3,i-time_horizon:i+time_horizon)).^2;
            N = 2*time_horizon + 1;
        elseif i > length(x(1,:)) - time_horizon
            d = (x(3,i) - x_star(3,i-time_horizon:length(x(1,:)))).^2;
            N = time_horizon + (length(x(1,:)) - i) + 1;
        end
        
        % Define Q as a diagonal matrix
        Q = diag(d);
        
        % Define constraint matrix A and vector b
        A = ones(1, N);
        b = 1;
        
        % Define binary constraints
        lb = zeros(N, 1);
        ub = ones(N, 1);
        
        % Solve MIQP with Gurobi
        model.Q = sparse(Q);
        model.A = sparse(A);
        model.rhs = b;
        model.sense = '=';
        model.lb = lb;
        model.ub = ub;
        model.vtype = 'B'; % Binary variables
        model.modelsense = 'min';
        
        % params.outputflag = 1;
        params.outputflag = 0;
        params.TimeLimit = 300;
        params.Threads = 2;

        result = gurobi(model, params);
        theta_error = [theta_error, sqrt(result.objval)]; % Path error
        fprintf('Set 3, Iteration %d from %d and solver time is %2.4f\n', i, length(x(1,:)), result.runtime);
    end

end