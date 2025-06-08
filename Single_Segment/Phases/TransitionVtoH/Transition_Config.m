function params = Transition_Config()
    %% Define parameters

    % Base parameters
    rho = 1.223;
    mass = 16.6;
    g = 9.81;
    S_ref = 0.5008;
    CL_max = 2.045; 
    Safety_Factor = 1.2;
    CL_max_ope = CL_max/Safety_Factor^2;

    % Prop system
    diameter = 0.8128;              % Propeller diameter [m]
    num_engines = 2;               % Number of engines
    rps_max = 78.125;        % Max rotations per second
    T_max_eng = 238.383;     % Max thrust per engine [N]
    P_max_eng = 6.72e3;     % Max power per engine [N]
    n_hov = 45.53;

    % Transition-specific parameters
    V_min_ope = sqrt(2*mass*g/(rho*S_ref*CL_max_ope));             % [m/s]
    epsilon_slope = deg2rad(2);  % Initial tilt rate [rad/s]
    alpha_opt_eff = deg2rad(2.66); %2.66
    phi_0 = deg2rad(90);
    alpha_mp = 15.9257*pi/180;% 16.4*pi/180;     %[rad] 16.4
    alpha_mi = -5*pi/180;       %[rad] -5
    tau_tr   = 12;              %[s]
    alpha_t  = 15.9257*pi/180; %12.9217*pi/180;  %[rad]
    V_f      = 19.35;              %[m/s]
    alpha_max_ope = deg2rad(15.9257);
    eta_e = 0.88*0.98*0.95;
    mass_batt= 3;                             % [kg]
    E_batt = 720e3*mass_batt;             % [J]
    
    % Solver settings
    tspan = [0 12];
    RelTol = 1e-6;

    % Plot options for transition corridor
    flag_corridor = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Define structures from parameters
    % Base parameters
    params.rho = rho;
    params.mass = mass;
    params.g = g;
    params.S_ref = S_ref;
    params.CL_max_ope = CL_max_ope;
    params.alpha_max_ope = alpha_max_ope;
    
    % Propulsion subsystem parameters
    params.prop.diameter = diameter;              
    params.prop.num_engines = num_engines;        
    params.prop.rps_max = rps_max;     
    params.prop.T_max_eng = T_max_eng;    
    params.prop.P_max_eng = P_max_eng;    
    params.prop.n_hov = n_hov;
    
    % Transition-specific parameters
    params.V_min_ope = V_min_ope;               
    params.epsilon_slope = epsilon_slope; 
    params.alpha_opt_eff = alpha_opt_eff; 
    params.phi_0 = phi_0;
    params.eta_e = eta_e;
    params.mass_batt= mass_batt;             
    params.E_batt = E_batt;
    params.alpha_mi = alpha_mi;
    params.alpha_t  = alpha_t;
    params.V_f      = V_f;
    params.alpha_mp = alpha_mp;
    params.tau_tr   = tau_tr; 

    
    % Solver settings
    params.tspan = tspan;
    params.RelTol = RelTol;

    % Plot options for transition corridor

    params.flag_corridor = flag_corridor;
end