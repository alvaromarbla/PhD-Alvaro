function [dy] = Phase2_DAE(t, y, params)
    % State variables: [x; V; E; α]
    % Algebraic equation: L + T·sinφ = W
    
    % Current controls
    epsilon = params.eps_final;
    alpha = y(4);
    
    % Aerodynamic forces
    L = CL_Model(alpha) * 0.5*params.rho*y(2)^2*params.S_ref;
    D = CD_Model(alpha) * 0.5*params.rho*y(2)^2*params.S_ref;
    
    % Propulsion forces (fixed at max RPM)
    [T, ~] = CT_Model(y(2), params.prop.rps_max, alpha, epsilon, params);
    P = CP_Model(y(2), params.prop.rps_max, alpha, epsilon, params);
    
    % Differential equations
    dy = zeros(4,1);
    dy(1) = y(2);                         % dx/dt
    dy(2) = (T*cos(epsilon + alpha) - D)/params.mass; % dV/dt
    dy(3) = P / params.eta_e;             % Energy consumption
    dy(4) = L + T*sin(epsilon + alpha) - params.mass*params.g;    % Placeholder for algebraic var
    

end