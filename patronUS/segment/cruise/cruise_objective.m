function J = cruise_objective(x, params)
    % Unpack variables: [V, \gamma, \alpha, \epsilon, n]
    V = x(1); alpha = x(3); epsilon = x(4); n = x(5); 

    % Power consumption
    
    CP = CP_Model(V, n, alpha, epsilon, params);
    P = params.prop.num_engines * CP * params.rho * n^3 * params.prop.diameter^5;
    
    J = -params.prop.eff*params.E_batt*V/(P);  % Minimize the minus range 
end