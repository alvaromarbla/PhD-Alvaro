function energy = Transition_Objective(epsilon_slope, base_params)
    try
        % Update tilt rate parameter
        params = base_params;
        params.epsilon_slope = epsilon_slope;
        
        % Run full transition simulation
        [~, energy] = Transition_Solver(params);
        
    catch ME
        warning('Failed for ε_slope=%.2f deg/s: %s',...
                rad2deg(epsilon_slope), ME.message);
        energy = inf;
    end
end