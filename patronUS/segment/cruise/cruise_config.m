function [params, bounds] = cruise_config (params)


%% This section only accounts for the PHASE SPECIFIC parameters.

% Cruise Segment.

% This accounts for segment specific params and bounds

mass_batt= 3;        % [kg]
tau = 0.2;  % Security factor for battery discharge
e0  = 720e3;  % Battery Energy density  [J/kg]
E_batt = e0*mass_batt*(1-tau); % [J]

params.E_batt = E_batt; % Cruise Energy


% Operational bounds (all values in SI units)
bounds.V.min = 10;                % Minimum airspeed [m/s]
bounds.V.max = 40;                % Maximum airspeed [m/s]

bounds.gamma.min = deg2rad(0);  % Minimum flight path angle [rad]
bounds.gamma.max = deg2rad(0);  % Maximum flight path angle [rad]

bounds.alpha.min = deg2rad(-30);  % Minimum angle of attack [rad]
bounds.alpha.max = deg2rad(45);   % Maximum angle of attack [rad]

bounds.epsilon.min = deg2rad(-90);% Minimum thrust vector angle [rad]
bounds.epsilon.max = deg2rad(90); % Maximum thrust vector angle [rad]

bounds.n.min = 20;                 % Minimum RPM
bounds.n.max = params.prop.max_rps; % Maximum RPM from physical limits





end 