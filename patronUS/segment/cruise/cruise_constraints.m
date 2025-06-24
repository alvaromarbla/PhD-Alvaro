function [c, ceq] = cruise_constraints(x, params,bounds)
    % Unpack variables
    V = x(1); alpha = x(3); epsilon = x(4); n = x(5); %
    
    % Reuse models from descent

    L = CL_wholeAC_1wing(alpha) * 0.5*params.rho*V^2*params.wing_area;
    D = CD_wholeAC_1wing(alpha) * 0.5*params.rho*V^2*params.wing_area;

    phi = alpha + epsilon;
    T = prop_CT_1365(V, n , phi, params);
    P = prop_CP_1365(V, n , phi, params);
    
    % Dynamics constraints (steady cruis)

    ceq = [T*cos(phi) - D ;   % Axial force balance
        L + T*sin(phi) - params.mass*params.g];  % Normal force balance

    % Bounds as inequality constraints
    c = [-V/(n*params.prop.diameter)  ; -T; -P ; T - params.prop.T_max_eng*params.prop.num_engines ;  P - params.prop.P_max_eng*params.prop.num_engines];      
end