function [c, ceq] = Climb_Constraints(x, params,bounds)
    % Unpack variables
    V = x(1); gamma = x(2); alpha = x(3); epsilon = x(4); n = params.n_ope;
    
    % Reuse models from descent
    [T, phi] = CT_Model(V, n, alpha, epsilon, params);
    CP = CP_Model(V, n, alpha, epsilon, params);
    L = CL_Model(alpha) * 0.5*params.rho*V^2*params.wing_area;
    D = CD_Model(alpha) * 0.5*params.rho*V^2*params.wing_area;
    P = params.prop.num_engines * CP * params.rho * n^3 * params.prop.diameter^5;
    t_climb = params.deltaH / (V*sin(gamma));
    % Dynamics constraints (steady climb)
    ceq = [T*cos(phi) - D - params.mass*params.g*sin(gamma);   % Axial force balance
        L + T*sin(phi) - params.mass*params.g*cos(gamma)];  % Normal force balance

    % Bounds as inequality constraints
    c = [V/(n*params.prop.diameter) - 1 ; -T; -P ; T - params.prop.T_max_eng*params.prop.num_engines ;  P - params.prop.P_max_eng*params.prop.num_engines; t_climb-600; -gamma];      
end