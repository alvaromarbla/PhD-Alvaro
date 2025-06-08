function J = Climb_Objective(x, params)
    % Unpack variables: [V, γ, α, ε, n]
    V = x(1); gamma = x(2); alpha = x(3); epsilon = x(4); n = params.n_ope;
    
    % Power consumption
    CP = CP_Model(V, n, alpha, epsilon, params);
    
    % Energy calculation
    t_climb = params.deltaH / (V*sin(gamma));
    Energy = params.prop.num_engines * CP * params.rho * n^3 * params.prop.diameter^5 * t_climb/params.prop.eff;
    
    % Horizontal distance benefit
    x_benefit = params.deltaH * cot(gamma) * params.C_cost;
    
    J = Energy- x_benefit;  % Minimize this
end