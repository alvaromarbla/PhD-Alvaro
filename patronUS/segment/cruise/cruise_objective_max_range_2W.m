function obj_fun = cruise_objective_max_range_2W(X_struct, models, params)
    % Unpack variables to match the 5-element vector: [V, gamma, alpha, epsilon, n]
    V       = X_struct.V;
    epsilon1 = X_struct.epsilon1;
    epsilon2 = X_struct.epsilon2;
    n1      = X_struct.n1;
    n2      = X_struct.n2;
    alpha   =  X_struct.alpha;
    %CP_lookup = models.CP_lookup;
    CT_lookup = models.CT_lookup; %%% Change
    % Power consumption models
    J1   = V / (n1 * params.prop.diameter);
    phi1 = alpha + epsilon1;
    J2   = V / (n2 * params.prop.diameter);
    phi2 = alpha + epsilon2;
    % CP1 = CP_lookup(J1,phi1);
    % CP2 = CP_lookup(J2,phi2); CHANGE THIS
    CT1 = CT_lookup(J1,phi1);
    CT2 = CT_lookup(J2,phi2);

    % P1 = params.prop.num_engines * CP1 * params.rho * (n1^3) * (params.prop.diameter^5);
    % P2 = params.prop.num_engines * CP2 * params.rho * (n2^3) * (params.prop.diameter^5);

    T1 = params.prop.num_engines * CT1 * params.rho * (n1^2) * (params.prop.diameter^4);
    T2 = params.prop.num_engines * CT2 * params.rho * (n2^2) * (params.prop.diameter^4);

    P1 = T1*V;
    P2 = T2*V;

    % Objective expression to minimize
    %obj_fun = -params.prop.eff * params.E_batt * V / (P);
    obj_fun = (P1+P2)/(V) * params.prop.diameter* params.prop.max_rps / params.prop.P_max_eng; % to make it non-dim
end