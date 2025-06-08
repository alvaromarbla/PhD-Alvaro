function [params, bounds] = Descent_Config()
    mass_batt= 3;                             % [kg]
    E_batt = 720e3*mass_batt;
    % Physical parameters
    params.mass = 16.6;               % [kg]
    params.wing_area = 0.5008;        % [m²]
    params.prop.diameter = 0.8128;    % [m]
    params.prop.num_engines = 2;
    params.prop.max_rps = 78.125;     % [rps] 
    params.prop.eff = 0.733;
    params.rho = 1.2133;              % [kg/m³]
    params.g = 9.81;                  % [m/s²]
    params.deltaH = 350;              % [m] Descent altitude
    params.alpha_mp = 16.4*pi/180;     %[rad] 16.4
    params.alpha_mi = -5*pi/180;       %[rad] -5
    % Phase-specific parameters
    params.phase = 'descent';
    params.C_cost = 20.039;
    params.mass_batt= mass_batt;             
    params.E_batt = E_batt; 
    
    % Operational bounds (all values in SI units)
    bounds.V.min = 10;                % Minimum airspeed [m/s]
    bounds.V.max = 40;                % Maximum airspeed [m/s]
    
    bounds.gamma.min = deg2rad(-90);             % Minimum flight path angle [rad]
    bounds.gamma.max = deg2rad(-1);   % Maximum flight path angle [rad]
    
    bounds.alpha.min = deg2rad(params.alpha_mi);  % Minimum angle of attack [rad]
    bounds.alpha.max = deg2rad(params.alpha_mp);   % Maximum angle of attack [rad]
    
    bounds.epsilon.min = deg2rad(-90);% Minimum thrust vector angle [rad]
    bounds.epsilon.max = deg2rad(90); % Maximum thrust vector angle [rad]
    
    bounds.n.min = 20;                 % Minimum RPM
    bounds.n.max = params.prop.max_rps;%params.prop.max_rps;%params.prop.max_rps; % Maximum RPM from physical limits
end