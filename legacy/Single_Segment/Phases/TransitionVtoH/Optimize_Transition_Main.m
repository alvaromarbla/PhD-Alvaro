function [opt_slope, min_energy] = Optimize_Transition_Main()
    % Initialize base parameters
    base_params = Transition_Config();  % Contains aircraft configuration
    
    % Define optimization bounds [deg/s]
    lower_bound = deg2rad(0.5);  % 0.5 deg/s
    upper_bound = deg2rad(60);   % 60 deg/s
    
    % Configure optimization
    options = optimset('Display', 'iter',...
                      'TolX', 1e-3,...
                      'OutputFcn', @optim_logger);
    
    % Run optimization
    [opt_slope, min_energy] = fminbnd(@(es) Transition_Objective(es, base_params),...
                                      lower_bound, upper_bound, options);
    
    % Final analysis
    fprintf('\nOptimal Tilt Rate: %.2f deg/s\n', rad2deg(opt_slope));
    fprintf('Minimum Energy: %.2f kJ\n', min_energy/1000);
    Plot_Optimal_Transition(opt_slope, base_params);
end

function stop = optim_logger(x, optimValues, state)
    % Custom output function to monitor optimization progress
    persistent history;
    
    stop = false; % Never stop early
    if strcmp(state, 'init')
        history = [];
        fprintf('%-12s %-12s %-12s\n', 'Iteration', 'Slope(deg/s)', 'Energy(kJ)');
    end
    
    if strcmp(state, 'iter')
        % Store and display current iteration
        history = [history; x optimValues.fval];
        
        fprintf('%-12d %-12.2f %-12.1f\n',...
                optimValues.iteration,...
                rad2deg(x),...
                optimValues.fval/1000);
    end
end