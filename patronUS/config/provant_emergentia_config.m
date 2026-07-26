function params = provant_emergentia_config

    mass_batt= 3;                             % [kg]
    E_batt = 720e3*mass_batt;

    % Physical parameters
    params.mass = 16.6;               % [kg]
    params.wing_area = 0.5008;        % [m²]
    params.prop.diameter = 0.8128;    % [m]
    params.prop.num_engines = 2;
    params.prop.max_rps = 78.125;
    params.prop.T_max_eng = 238.383;  % [N] Max thrust per engine
    params.prop.P_max_eng = 6.7e3;    % [W] Max power per engine
    params.prop.eff = 0.733;        % [-] Efficiency of engines
    params.rho = 1.2133;              % [kg/m³]
    params.g = 9.81;                  % [m/s²]
    params.mass_batt= mass_batt;             
    params.E_batt = E_batt; 
end
  
