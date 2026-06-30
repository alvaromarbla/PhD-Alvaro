function [sol, fval] = Main_Descent_Optimizer()
    % Initialize configuration
    [params, bounds] = Descent_Config();
    
    % Optimization variables: [V, gamma, alpha, epsilon, n]
    vars0 = [19, deg2rad(-5), deg2rad(15), deg2rad(35), 25]; % Initial guess
    

    % Create bound vectors using configured constraints
    lb = [bounds.V.min, bounds.gamma.min, bounds.alpha.min, ...
          bounds.epsilon.min, bounds.n.min];
    ub = [bounds.V.max, bounds.gamma.max, bounds.alpha.max, ...
          bounds.epsilon.max, bounds.n.max];
    
    % Configure solver
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'iter', ...
        'MaxFunctionEvaluations', 1e5,'ConstraintTolerance',1e-10);
    
    % Solve
    [sol, fval] = fmincon(...
        @(x) MinusCotGamma(x), ... % J_descent(x, params) // MinusCotGamma(x)
        vars0, [], [], [], [], lb, ub, ...
        @(x) Dynamics_Constraints(x, params), ...
        options);
    
    % Post-process
    Visualization(sol, params);
    %Save_Results(sol, params);
end