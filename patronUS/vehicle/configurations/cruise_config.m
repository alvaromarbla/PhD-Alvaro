function [params, bounds] = cruise_config(AC_version)

import casadi.*

if AC_version == "1WING"
    mass_batt= 3;                                  % [kg]
    E_batt = 720e3*mass_batt;

    % Physical parameters
    params.mass = 16.6;                            % [kg]
    params.wing_area = 0.5008;                     % [m²]
    params.nacelle_area = 0.3;                     % [m^2] to be decided
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
    bounds.V.min = 10;                             % Minimum airspeed [m/s]
    bounds.V.max = 40;                             % Maximum airspeed [m/s]

    bounds.gamma.min = deg2rad(-0.5);   % Minimum flight path angle [rad]
    bounds.gamma.max = deg2rad(0.5);    % Maximum flight path angle [rad]

    bounds.alpha.min = deg2rad(-30);    % Minimum angle of attack [rad]
    bounds.alpha.max = deg2rad(45);     % Maximum angle of attack [rad]

    bounds.epsilon.min = deg2rad(0);% Minimum thrust vector angle [rad]
    bounds.epsilon.max = deg2rad(90); % Maximum thrust vector angle [rad]

    bounds.n.min = 25;                              % Minimum RPM
    bounds.n.max = params.prop.max_rps; % Maximum RPM from physical limits

elseif AC_version == "2WING"

    mass_batt= 3;                                  % [kg]
    E_batt = 720e3*mass_batt;

    % Physical parameters
    params.mass = 16.6;                            % [kg]
    params.wing_area = 0.5008;                     % [m²]
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
    bounds.V.min = 10;                             % Minimum airspeed [m/s]
    bounds.V.max = 40;                             % Maximum airspeed [m/s]

    bounds.gamma.min = deg2rad(-0.5);   % Minimum flight path angle [rad]
    bounds.gamma.max = deg2rad(0.5);    % Maximum flight path angle [rad]

    bounds.alpha.min = deg2rad(-30);    % Minimum angle of attack [rad]
    bounds.alpha.max = deg2rad(45);     % Maximum angle of attack [rad]

    bounds.epsilon.min = deg2rad(0);% Minimum thrust vector angle [rad]
    bounds.epsilon.max = deg2rad(90); % Maximum thrust vector angle [rad]

    bounds.n.min = 25;                              % Minimum RPM
    bounds.n.max = params.prop.max_rps; % Maximum RPM from physical limits

else
    warning('AC version unknown')

end