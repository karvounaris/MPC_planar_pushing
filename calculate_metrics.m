%============================================================================%
% Function to calculate performance, control effort, robustness, and computation metrics
%
% Inputs:
% - x: Matrix of states (size: 4 x time steps)
% - x_star: Matrix of reference states (size: 4 x time steps)
% - u: Matrix of control inputs (size: 3 x time steps)
% - time: Vector of time steps
% - solver_times: Vector of solver computation times (size: 1 x time steps)
%
% Outputs:
% - metrics: Struct containing the calculated metrics
%============================================================================%
function metrics = calculate_metrics(x, x_star, u, solver_times, N)
    % Performance Evaluation Metrics
    % Errors for the first two states
    state_errors = x(1:2, end) - x_star(1:2, end-N);
    
    % RMSE: Aggregate across the first two states
    RMSE = sqrt(mean(state_errors.^2)); % Single scalar for combined error
    
    % Maximum deviation: Aggregate across the first two states
    max_deviation = max(abs(state_errors)); % Single scalar for combined deviation

    % Maximum magnitude of each control input component
    max_control_magnitude = max(abs(u), [], 2);

    % Terminal state deviation: Norm of the terminal state error
    terminal_state_deviation = norm(state_errors); % Single scalar as the Euclidean norm

    % Computation Metrics
    % Average solver time
    avg_solver_time = mean(solver_times);

    % Pack results into a struct
    metrics.RMSE = RMSE; % Scalar
    metrics.MaxDeviation = max_deviation; % Scalar
    metrics.MaxControlMagnitude = max_control_magnitude;
    metrics.TerminalStateDeviation = terminal_state_deviation; % Scalar
    metrics.AvgSolverTime = avg_solver_time;
end

