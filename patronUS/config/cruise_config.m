function [params, bounds] = cruise_config(AC_version)

import casadi.*

if AC_version == "1Wing"
    mass_batt= 3;                                  % [kg]
    E_batt = 720e3*mass_batt;

    % Physical parameters
    params.mass = 16.6;                            % [kg]
    params.wing_area = 0.5008;                     % [m²]
    params.nacelle_area = 0.1;                     % [m^2] to be decided
    params.prop.diameter = 0.8128*1;                % [m]
    params.prop.num_engines = 2;
    params.prop.max_rps = 78.125;
    params.prop.T_max_eng = 238.383;               % [N] Max thrust per engine
    params.prop.P_max_eng = 6.7e3;                 % [W] Max power per engine
    params.prop.eff = 0.8192;%0.733;                    % [-] Efficiency of engines
    params.rho = 1.2133;                           % [kg/m³]
    params.g = 9.81;                               % [m/s²]
    params.deltaH = 350;                           % [m] Climb altitude
    params.mass_batt= mass_batt;
    params.E_batt = E_batt;

    % Phase-specific parameters
    params.phase = 'cruise';
    params.deltaH = 350;                           % [m] Target altitude gain
    params.delta0 = 350;                            % [m] Initial altitude

    % Operational bounds (all values in SI units)
    bounds.vars.V       = struct('min', 10,                'max', 40,               'init', 24.225);
    bounds.vars.gamma   = struct('min', deg2rad(-0.5),'max', deg2rad(0.5),           'init', deg2rad(0));
    bounds.vars.alpha   = struct('min', deg2rad(-30), 'max', deg2rad(45),            'init', 0.0464);
    bounds.vars.epsilon = struct('min', deg2rad(0),   'max', deg2rad(90),            'init', deg2rad(60));
    bounds.vars.n       = struct('min', 25,                 'max', params.prop.max_rps, 'init', 25.214);

elseif AC_version == "2Wings"

    mass_batt= 3;                                  % [kg]
    E_batt = 720e3*mass_batt;

    % Physical parameters
    params.mass = 255.4;                           % [kg]
    params.Iyy  = 98.19;                            % [kg m^2]
    params.Iyy_nac = 0.06469;                       % [kg m^2]

    params.wing_area = 0.48186;                     % [m²]
    params.tailwing_area = 0.48186;                 % [m²]
    params.fus_area  = 0.4094;                      % [m²]
    params.fus_Lref  = 1.75;                        % [m]

    params.prop.diameter = 0.718*1;                 % [m]
    params.prop.num_engines = 2;
    params.prop.max_rps = 78.125;
    params.prop.T_max_eng = 238.383;               % [N] Max thrust per engine
    params.prop.P_max_eng = 6.7e3;                 % [W] Max power per engine
    params.prop.eff = 0.8192;%0.733;                    % [-] Efficiency of engines
    params.prop.nacelle_area = 0.01774;             % [m²]
    params.rho = 1.2133;                           % [kg/m³]
    params.g = 9.81;                               % [m/s²]
    params.deltaH = 350;                           % [m] Climb altitude
    params.mass_batt = mass_batt;
    params.mass_eng  = 14.35;
    params.E_batt = E_batt;

    params.geo.c_wing = 0.23;                       % [m] Mean wing chord
    params.geo.xw     = 0.707;                      % [m] X distance
    params.geo.zw     = 0.0805;
    params.geo.xfus   = 0.707/2; %%%%% THESE ARE TO CHANGE
    params.geo.zfus   = 0.0805/2; %%%%% THESE ARE TO CHANGE
    params.geo.xtw    = 0.619;
    params.geo.ztw    = 0.1935;

    % Phase-specific parameters
    params.phase = 'cruise';
    params.deltaH = 350;                           % [m] Target altitude gain
    params.delta0 = 350;                            % [m] Initial altitude

    % Operational bounds (all values in SI units)
    bounds.vars.V        = struct('min', 10,               'max', 40,                    'init', 30.0);
    bounds.vars.gamma    = struct('min', deg2rad(0),'max', deg2rad(0),               'init', deg2rad(0));
    bounds.vars.theta    = struct('min', deg2rad(-30), 'max', deg2rad(45),                'init', deg2rad(0.1));
    bounds.vars.alpha    = struct('min', deg2rad(-30), 'max', deg2rad(45),                'init', 0.1);
    bounds.vars.epsilon1 = struct('min', deg2rad(0),   'max', deg2rad(90),                'init', deg2rad(4));
    bounds.vars.epsilon2 = struct('min', deg2rad(0),   'max', deg2rad(90),                'init', deg2rad(4));
    bounds.vars.q        = struct('min', -0.6,             'max', 0.6,                    'init', 0); % Rad/s

    % Controls
    bounds.vars.n1   = struct('min', 25,   'max', params.prop.max_rps, 'init', 40.0);
    bounds.vars.n2   = struct('min', 25,   'max', params.prop.max_rps, 'init', 40.0);
    % bounds.vars.tau1 = struct('min', -5,   'max', 5,                   'init', 0.0);
    % bounds.vars.tau2 = struct('min', -5,   'max', 5,                   'init', 0.0);
    bounds.vars.deltae = struct('min', deg2rad(-20), 'max', deg2rad(20),                  'init', 0);

else
    warning('AC version unknown')

end