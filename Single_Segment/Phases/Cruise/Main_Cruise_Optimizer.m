function [sol, fval] = Main_Cruise_Optimizer

    % Initialize configuration
    [params, bounds] = Cruise_Config;

    % Optimization variables: [V, gamma, alpha, epsilon, n]
    %vars0 = [25, deg2rad(36), deg2rad(-5), deg2rad(0), 55]; % Initial guess
    vars0 = [24.225109, 0.046444, deg2rad(60), 25.214554]; % Initial guess 1.0593

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


    % Solve
    [sol, fval] = fmincon(...
        @(x) Cruise_Objective(x, params), ...
        vars0, [], [], [], [], lb, ub, ...
        @(x) Cruise_Constraints(x, params,bounds), ...
        options);

    V_sol = sol(1);
    n_sol = sol(4);
    alpha_sol = sol(2);
    epsilon_sol = sol(3);
    [T_sol, phi_sol] = CT_Model(V_sol, n_sol, alpha_sol, epsilon_sol, params)
    CP = CP_Model(V_sol, n_sol, alpha_sol, epsilon_sol, params);
    CL_Model(alpha_sol)
    CD_Model(alpha_sol)
    CL_Model(alpha_sol)/ CD_Model(alpha_sol)
    L_sol = CL_Model(alpha_sol) * 0.5*params.rho*V_sol^2*params.wing_area
    D_sol = CD_Model(alpha_sol) * 0.5*params.rho*V_sol^2*params.wing_area
    P_sol = params.prop.num_engines * CP * params.rho * n_sol^3 * params.prop.diameter^5
    % Post-process
    %Cruise_2D_plot(sol,fval,params);
    %Jsol = Cruise_Objective(sol,params)
end