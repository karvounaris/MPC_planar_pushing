%========================================================================%
% Function to perform a grid search for MPC parameters and evaluate metrics
%
% Inputs:
% - Q_values: Cell array of candidate Q matrices
% - QN_values: Cell array of candidate QN matrices
% - R_values: Cell array of candidate R matrices
% - mpc_timesteps: Array of candidate MPC timesteps
% - N_values: Array of candidate N values
% - duration: Simulation duration
% - timestep: Simulation timestep
%
% Outputs:
% - results: Struct array containing parameter sets and their corresponding metrics
%========================================================================%

function results = grid_search(x_star, u_star, Q_values, QN_values, R_values, mpc_timesteps, ...
                               N_values, duration, timestep, mu, L, radius, len, mass, x_initial, I_object)

    % Initialize results storage
    results = [];
    idx = 1; % Index for storing results

    % Iterate over all combinations of Q, QN, R, and mpc_timestep
    for qi = 1:length(Q_values)
        for qni = 1:length(QN_values)
            for ri = 1:length(R_values)
                for mpc_idx = 1:length(mpc_timesteps)
                    for Ni = 1:length(N_values)
                        % Set current parameters
                        Q = Q_values{qi};
                        QN = QN_values{qni};
                        R = R_values{ri};
                        mpc_timestep = mpc_timesteps(mpc_idx);
                        N = N_values(Ni);
    
                        % Display progress
                        fprintf('Running simulation with Q[%d], QN[%d], R[%d], mpc_timestep[%d] and N[%d]\n', ...
                            qi, qni, ri, mpc_idx, Ni);
    
                        try
                            % Call the simulation function
                            metrics = run_simulation(Q, QN, R, mpc_timestep, duration, timestep, ...
                                                     x_star, u_star, mu, L, radius, len, N, mass, x_initial, I_object);

                            % Store results
                            results(idx).Q = Q;
                            results(idx).QN = QN;
                            results(idx).R = R;
                            results(idx).mpc_timestep = mpc_timestep;
                            results(idx).N = N;
                            results(idx).RMSE = metrics.RMSE;
                            results(idx).TerminalStateDeviation = metrics.TerminalStateDeviation;
                            results(idx).AvgSolverTime = metrics.AvgSolverTime;
                            results(idx).MaxDeviation = metrics.MaxDeviation;
                            results(idx).MaxControlMagnitude = metrics.MaxControlMagnitude;
                            idx = idx + 1;
                        catch ME
                            fprintf('Error in simulation (Q[%d], QN[%d], R[%d], timestep[%d], N[%d]): %s\n', ...
                                qi, qni, ri, mpc_idx, Ni, ME.message);
                        end
                    end
                end
            end
        end
    end
end
