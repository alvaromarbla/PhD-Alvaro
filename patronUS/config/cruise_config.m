function [params, bounds,models] = cruise_config(versions)

import casadi.*

models = models_config(versions);


if versions.AC_version == "1Wing"
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

elseif versions.AC_version == "2Wings"

    mass_batt= 3;                                  % [kg]
    E_batt = 720e3*mass_batt;

    % Physical parameters
    params.mass = 25.54;         %  255.4 says CAD?                 % [kg]
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
    bounds.vars.V        = struct('min', 10,               'max', 40,                    'init', 22.0);
    bounds.vars.gamma    = struct('min', deg2rad(0),'max', deg2rad(0),               'init', deg2rad(0));
    bounds.vars.theta    = struct('min', deg2rad(-30), 'max', deg2rad(45),                'init', deg2rad(14));
    %bounds.vars.alpha    = struct('min', deg2rad(-30), 'max', deg2rad(45),                'init', deg2rad(3));
    bounds.vars.epsilon1 = struct('min', deg2rad(0),   'max', deg2rad(90),                'init', deg2rad(40));
    bounds.vars.epsilon2 = struct('min', deg2rad(0),   'max', deg2rad(90),                'init', deg2rad(40));
    bounds.vars.q        = struct('min', -0.6,             'max', 0.6,                    'init', 0); % Rad/s

    % Controls
    bounds.vars.n1   = struct('min', 18,   'max', params.prop.max_rps, 'init', 25.0);
    bounds.vars.n2   = struct('min', 18,   'max', params.prop.max_rps, 'init', 25.0);
    % bounds.vars.tau1 = struct('min', -5,   'max', 5,                   'init', 0.0);
    % bounds.vars.tau2 = struct('min', -5,   'max', 5,                   'init', 0.0);
    bounds.vars.deltae = struct('min', deg2rad(-20), 'max', deg2rad(20),                  'init', 0);

else
    warning('AC version unknown')

end

%% Battery model: ideal energy store or Shepherd pack (versions.Batt_version)
% The Shepherd pack replaces the reference battery above (mass and energy)

params = battery_config(params, versions);

%% Time-dependent cruise (Shepherd battery, 2Wings)
% With the Shepherd model the range for a given charge is no longer
% eta*E*V/P: the current needed for a power grows as SOC drops, so the
% cruise is an optimal control problem in time (main_cruise_optimizer):
% the trim of bounds.vars is imposed at every node (quasi-steady, the
% variables become controls), and the states are the range x, the bus
% energy E and SOC, from SOC0 down to SOC_f, with tf free.

if strcmp(params.batt.model, 'Shepherd') && versions.AC_version == "2Wings"

    params.cruise.SOC_f   = 0.3;                   % [-] End of cruise: 10 % descent, transition and landing + 20 % safety reserve
    params.ocp.N          = 20;                    % Number of collocation intervals
    params.ocp.reg_weight = 1e-3;                  % [-] Control-rate regulariser (control_rate_penalty), 0 = off
    params.ocp.x_ref      = 100e3;                 % [m] Range scaling of the objective

    SOC0 = params.batt.SOC0;
    bounds.states.x   = struct('min', 0,    'max', 500e3,         'x0', 0,    'xf', NaN,                 'init', [0 60e3]);  % Range
    bounds.states.E   = struct('min', 0,    'max', params.E_batt, 'x0', 0,    'xf', NaN,                 'init', [0 0.7*params.E_batt]); % Bus energy drawn
    bounds.states.SOC = struct('min', 0.05, 'max', 1,             'x0', SOC0, 'xf', params.cruise.SOC_f, 'init', [SOC0 params.cruise.SOC_f]);

    bounds.controls = bounds.vars;                 % Trim variables, free at every node

    % scale: tf/scale is the decision variable (cruise lasts ~1 h)
    bounds.time.tf = struct('min', 60, 'max', 6*3600, 'init', 3600, 'scale', 3600);

end