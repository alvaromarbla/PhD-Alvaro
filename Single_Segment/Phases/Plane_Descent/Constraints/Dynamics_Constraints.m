function [c, ceq] = Dynamics_Constraints(x, params)
    % Extract variables
    V = x(1);
    gamma = x(2);
    alpha = x(3);
    epsilon = x(4);
    n = x(5);
    
    % Power constraint (C_P = 0)
    CP = CP_Model(V, n, alpha, epsilon, params);
    %CP0 = CP_Model(0, n, 0, 0, params);

    % Force calculations
    L = CL_Model(alpha) * 0.5 * params.rho * V^2 * params.wing_area;
    D = CD_Model(alpha) * 0.5 * params.rho * V^2 * params.wing_area;
    P = params.prop.num_engines * params.rho * n^3 * params.prop.diameter^5 *CP;
    t_climb = params.deltaH / (V*sin(gamma));
    Energy = params.prop.num_engines * CP * params.rho * n^3 * params.prop.diameter^5 * t_climb/params.prop.eff;
    [T, phi] = CT_Model(V, n, alpha, epsilon, params);
    
    
    % Dynamic constraints from problem formulation
    ceq = [T*cos(phi) - D - params.mass*params.g*sin(gamma);
        L + T*sin(phi) - params.mass*params.g*cos(gamma); % Power consumption = 0
    T];
   % T*cos(phi) - D - params.mass*params.g*sin(gamma);  % V-dot = 0
        %L + T*sin(phi) - params.mass*params.g*cos(gamma);  % gamma-dot = 0
    % Additional inequality constraints
    c = [V/(n*params.prop.diameter)-1.2, -P, +alpha+epsilon -pi/2 , Energy-params.E_batt*(0.1)];%;

   % - D - params.mass*params.g*sin(gamma);
    %    L - params.mass*params.g*cos(gamma)
end