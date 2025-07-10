function [Xsol,Xfval] = main_cruise 
    
    % Initialize configuration
    % Plane config 

    params = provant_emergentia_config;

    % Phase specific config
    [params, bounds] = cruise_config(params);

    % Optimization variables: [V, gamma, alpha, epsilon, n]
    

    % Create bound vectors using configured constraints
    lb = [bounds.V.min, bounds.gamma.min ,bounds.alpha.min, ...
          bounds.epsilon.min, bounds.n.min]';
    ub = [bounds.V.max,bounds.gamma.min  ,bounds.alpha.max, ...
          bounds.epsilon.max, bounds.n.max]';

    % Initial guess
    N_sample = 50; % In case we want to do a Latin Hypercube Sampling
    initial_cond = cruise_initial_cond (lb, ub, N_sample);

    % Configure solver
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'iter', ...
        'MaxFunctionEvaluations', 1e4,'ConstraintTolerance',1e-9);

    Xsol = zeros(length(ub),N_sample);
    Xfval = zeros(1,N_sample);
    k_V = 1; % Index for succesful runs;
    for kk = 1: N_sample
    % Solve %% Extract more information!!!
    [sol, fval,exitflag,output,lambda,grad,hessian] = fmincon(...
        @(x) cruise_objective(x, params), ...
        initial_cond(:,kk), [], [], [], [], lb, ub, ...
        @(x) cruise_constraints(x, params,bounds), ...
        options);

    if exitflag == 1
        Xsol(:,k_V) = sol;
        Xfval(k_V) = fval; 
        k_V = k_V +1;
    end 
   
    end

    % Post-process

end