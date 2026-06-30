function [dy] = Phase1_DAE(t, y, params)
    % State variables: [x; V; E; n]
    % Algebraic equation: L + T·sinφ = W
    alpha_mi = params.alpha_mi;
    alpha_t  = params.alpha_t;
    V_f      = params.V_f;
    tau_tr   = params.tau_tr;
    % Current controls
    %if y(2)<V_f
    %     alp  = alpha_mi + (alpha_t-alpha_mi)*(3*(V/V_f)^2-2*(V/V_f)^3);
        %alpha  = alpha_mi*(1-exp(-4*t)) + (alpha_t-alpha_mi)*(exp(y(2)/4)/exp(V_f/4));
    alpha  = alpha_t + (alpha_mi-alpha_t)*exp(-5*t/tau_tr);
        %epsilon = 0.1*pi/2*exp(-4*t) + 0.9*pi/2*(1-(exp(y(2)/4)/exp(V_f/4)));
    epsilon = pi/2*(1-t/tau_tr).^2;
    %else
        % alpha  = alpha_t;
        % epsilon = 0;
    %end

    % alpha = params.alpha_opt_eff;
    % t_est = params.phi_0/params.epsilon_slope; %Estimated time for the engines to finish tilt from the start  
    % epsilon = params.phi_0 - params.epsilon_slope*t - alpha + (params.epsilon_slope*t - params.phi_0 )*heaviside(t-t_est);
    
    % Aerodynamic forces
    L = CL_Model(alpha) * 0.5*params.rho*y(2)^2*params.S_ref;
    D = CD_Model(alpha) * 0.5*params.rho*y(2)^2*params.S_ref;
    
    % Propulsion forces
    [T, ~] = CT_Model(y(2), y(4), alpha, epsilon, params);
    P = CP_Model(y(2), y(4), alpha, epsilon, params)*params.rho*params.prop.num_engines*y(4).^3*params.prop.diameter^5;
    
    % Differential equations
    dy = zeros(4,1);
    dy(1) = y(2);                                                           % dx/dt
    dy(2) = (T*cos(epsilon + alpha) - D)/params.mass;                       % dV/dt
    dy(3) = P / params.eta_e;                                               % Energy consumption
    dy(4) = (L + T*sin(epsilon + alpha))/params.mass - params.g;              % Placeholder for algebraic var
    
   % if y(4)> params.prop.rps_max
    %    dy(4) = params.prop.rps_max;
    %end
end