function obj_fun = cruise_objective_max_range_1W(X_struct, models, params)
    % Unpack variables to match the 5-element vector: [V, gamma, alpha, epsilon, n]
    V       = X_struct.V;
    alpha   = X_struct.alpha;
    epsilon = X_struct.epsilon;
    n       = X_struct.n;

    CP_lookup = models.CP_lookup;     
    % Power consumption model
    J   = V / (n * params.prop.diameter);
    phi = alpha + epsilon;
    CP = CP_lookup(J,phi);

  
    P = params.prop.num_engines * CP * params.rho * (n^3) * (params.prop.diameter^5);

    % Objective expression to minimize
    %obj_fun = -params.prop.eff * params.E_batt * V / (P);
    obj_fun = P/(V) * params.prop.diameter* params.prop.max_rps / params.prop.P_max_eng; % to make it non-dim
end