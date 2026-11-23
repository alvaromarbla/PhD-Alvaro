function params = battery_config(params, versions)
% Battery model of the aircraft, shared by all segments (called from
% cruise_config, so vtol_config inherits it).
%   versions.Batt_version: ' ' (or missing) = ideal energy store (state E
%   only, reference battery params.mass_batt / params.E_batt of
%   cruise_config); 'Shepherd' = Shepherd pack (battery_shepherd), which
%   adds the state of charge SOC as a state in the trajectory segments.

batt_version = '';
if isfield(versions, 'Batt_version'), batt_version = strtrim(char(versions.Batt_version)); end

switch batt_version
    case ''
        params.batt.model = 'ideal';
    case 'Shepherd'
        % Pack sized by main_battery_sizing: P42A-class 21700 cells, 18S6P
        % (108 cells, 7.56 kg, ~1.68 kWh), for the take-off peak power
        % +20 % above the cut-off and within 45 A per cell.
        % battery_cell('V2') gives the ShepherdModel_V2 energy cell.
        params.batt = battery_pack(battery_cell('P42A'), 18, 6);
        params.batt.SOC0     = 1;                % [-]  Full battery at the start of the segment
        params.batt.disc_eps = 1;                % [V^2] Smoothing of max(disc, 0) in the current
        % The pack replaces the reference battery of cruise_config in the aircraft mass
        params.mass          = params.mass - params.mass_batt + params.batt.mass;
        params.mass_batt     = params.batt.mass;
        params.E_batt        = params.batt.E0*params.batt.Q;   % [J] Nominal pack energy (scales E)
    otherwise
        error('battery_config:unknown', 'Unknown Batt_version "%s" (use '' '' or ''Shepherd'').', batt_version);
end

end
