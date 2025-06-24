function J = cruise_objective(x, params)
    % Unpack variables: [V, \gamma, \alpha, \epsilon, n]
    V = x(1); alpha = x(3); epsilon = x(4); n = x(5); 

    %% Power consumption
    
    phi = alpha + epsilon;
    P = prop_CP_1365(V, n , phi, params);

    %% Objective Function
    J = -params.prop.eff*params.E_batt*V/(P);  % Minimize the minus range 
end