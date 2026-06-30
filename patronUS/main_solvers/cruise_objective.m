function obj_fun = cruise_objective(x, CP_lookup, params)
    % Unpack variables to match the 5-element vector: [V, gamma, alpha, epsilon, n]
    V       = x(1);
    alpha   = x(3);
    epsilon = x(4);
    n       = x(5);
    % Power consumption model
    J   = V / (n * params.prop.diameter);
    phi = alpha + epsilon;
    CP = CP_lookup(J,phi);

    % Core power calculation (CasADi natively handles ^ and * operators)
    P = params.prop.num_engines * CP * params.rho * (n^3) * (params.prop.diameter^5);

    % Objective expression to minimize
    %obj_fun = -params.prop.eff * params.E_batt * V / (P);
    obj_fun = P/(V*params.prop.eff * params.E_batt);
end