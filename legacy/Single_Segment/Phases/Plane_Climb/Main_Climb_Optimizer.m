function [sol, fval] = Main_Climb_Optimizer(C_cost)
    % Initialize configuration
    [params, bounds] = Climb_Config(C_cost);
    
    % Optimization variables: [V, gamma, alpha, epsilon, n]
    %vars0 = [25, deg2rad(36), deg2rad(-5), deg2rad(0), 55]; % Initial guess
    vars0 = [22, deg2rad(13), deg2rad(0), deg2rad(56)]; % Initial guess

    % Create bound vectors using configured constraints
    lb = [bounds.V.min, bounds.gamma.min, bounds.alpha.min, ...
          bounds.epsilon.min];
    ub = [bounds.V.max, bounds.gamma.max, bounds.alpha.max, ...
          bounds.epsilon.max];
    
    % Configure solver
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'iter', ...
        'MaxFunctionEvaluations', 1e4,'ConstraintTolerance',1e-9);
    
    % Solve
    [sol, fval] = fmincon(...
        @(x) Climb_Objective(x, params), ...
        vars0, [], [], [], [], lb, ub, ...
        @(x) Climb_Constraints(x, params,bounds), ...
        options);
    
    % Post-process
   
    %Visualization(sol, params);
    %Save_Results(sol, params);
end