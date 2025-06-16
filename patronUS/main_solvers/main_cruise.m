function [sol, fval] = main_cruise %Solve cruise?
    
    
    % Initialize configuration
    % Plane config 

    params = provant_emergentia_config;
    [params, bounds] = cruise_config(params);

    % Optimization variables: [V, gamma, alpha, epsilon, n]
    
    % Initial guess
    previous_cond = [0,0,0,0,0];
    initial_cond = cruise_initial_cond (previous_cond, params);

    % Create bound vectors using configured constraints
    lb = [bounds.V.min, bounds.alpha.min, ...
          bounds.epsilon.min, bounds.n.min];
    ub = [bounds.V.max, bounds.alpha.max, ...
          bounds.epsilon.max, bounds.n.max];

    % Configure solver
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'iter', ...
        'MaxFunctionEvaluations', 1e4,'ConstraintTolerance',1e-9);


    % Solve %% Extract more information!!!
    [sol, fval] = fmincon(...
        @(x) Cruise_Objective(x, params), ...
        initial_cond, [], [], [], [], lb, ub, ...
        @(x) Cruise_Constraints(x, params,bounds), ...
        options);

    % Post-process

end