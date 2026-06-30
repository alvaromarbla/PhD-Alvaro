function [c, ceq] = Cruise_Constraints(x, params,bounds)
    % Unpack variables
    V = x(1); alpha = x(2); epsilon = x(3); n = x(4); %
    
    % Reuse models from descent
    [T, phi] = CT_Model(V, n, alpha, epsilon, params);
    CP = CP_Model(V, n, alpha, epsilon, params);
    L = CL_Model(alpha) * 0.5*params.rho*V^2*params.wing_area;
    D = CD_Model(alpha) * 0.5*params.rho*V^2*params.wing_area;
    P = params.prop.num_engines * CP * params.rho * n^3 * params.prop.diameter^5;
    % Dynamics constraints (steady cruis)

    ceq = [T*cos(phi) - D ;   % Axial force balance
        L + T*sin(phi) - params.mass*params.g];  % Normal force balance

    % Bounds as inequality constraints
    c = [-V/(n*params.prop.diameter)  ; -T; -P ; T - params.prop.T_max_eng*params.prop.num_engines ;  P - params.prop.P_max_eng*params.prop.num_engines];      
end