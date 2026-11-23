%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% PATRONUS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Battery pack sizing for the VTOL take-off
% By: Alvaro Martinez Blanco (2026)
% The pack must deliver the peak bus power of the minimum-energy take-off
% with an ideal battery, plus a margin, above its cut-off voltage, and
% within the continuous current rating of the cell (battery_size_pack).
% The pack replaces the reference battery (params.mass_batt), so its mass
% changes the take-off and its peak power: the take-off is solved again at
% each new mass until the number of parallel strings stops changing.
% The result goes into the Shepherd block of vtol_config.

close all
%% User input

cell_name = 'P42A';      % battery_cell
Ns        = 18;          % Cells in series
margin    = 0.2;         % Power margin above the take-off peak

versions.AC_version    = "2Wings";
versions.CL_version    = 'analytical';
versions.CD_version    = 'analytical';
versions.CM_version    = 'analytical';
versions.CT_version    = 'analytical';
versions.CP_version    = 'analytical';
versions.CH_version    = 'zeros';
versions.Batt_version  = ' ';                 % Ideal battery: gives the power demand
versions.dynamics      = "VTOL_2Wings";
versions.constraints   = "VTOL_2Wings";
versions.objective     = "min_energy_vtol";
versions.initial_guess = "flat_2Wings";
versions.transcription = "trapezoidal";

%% Configuration

import casadi.*
[params0, bounds, models] = vtol_config(versions);
registry = patronus_registry();
validate_config(versions, registry);
fns.dynamics    = registry.dynamics(char(versions.dynamics));
fns.constraints = registry.constraints(char(versions.constraints));
fns.objective   = registry.objective(char(versions.objective));
guess_fcn       = registry.initial_guess(char(versions.initial_guess));
method          = registry.transcription(char(versions.transcription));

cell       = battery_cell(cell_name);
m_airframe = params0.mass - params0.mass_batt;    % Aircraft without the reference battery

%% Sizing loop (mass <-> peak power)

mass = params0.mass;
Np   = NaN;
sizing = struct('mass', {}, 'P_req', {}, 'Np', {}, 'P_cut', {}, 'I_cell', {}, 'U', {}, 'limit', {});
fprintf('\n  it   mass [kg]   tf [s]   peak P_bus [kW]   Np   pack [kg]   P_cut [kW]   I_cell [A]   U at peak [V]   limit\n');
for it = 1:10
    params = params0;
    params.mass = mass;

    [guess, ~] = select_initial_guess(guess_fcn, params.ocp.tf_sweep, fns.objective, fns.dynamics, ...
                                      fns.constraints, bounds, models, params);
    [res, info] = ocp_ipopt_stage(params, bounds, models, fns, method, guess);
    if ~info.ok
        warning('Take-off not converged at mass %.2f kg.', mass);
    end

    P_req = max(res.dbg.P_bus);
    [Np_new, s] = battery_size_pack(cell, Ns, P_req, margin);
    sizing(end+1) = struct('mass', mass, 'P_req', P_req, 'Np', Np_new, 'P_cut', s.P_cut, ...
                           'I_cell', s.I_cell, 'U', s.U, 'limit', s.limit); %#ok<SAGROW>
    fprintf('  %2d   %8.2f   %6.3f   %14.2f   %3d   %8.2f   %9.2f   %9.1f   %12.1f   %s\n', it, mass, res.tf, ...
            P_req/1e3, Np_new, s.mass, s.P_cut/1e3, s.I_cell, s.U, s.limit);

    if Np_new == Np
        break
    end
    Np   = Np_new;
    mass = m_airframe + s.mass;
end

pack = battery_pack(cell, Ns, Np);
fprintf('\nSized pack: %s %dS%dP, %d cells, %.2f kg, %.0f Wh nominal; aircraft mass %.2f kg\n', cell_name, Ns, Np, ...
        Ns*Np, pack.mass, pack.E0*pack.Q/3600, m_airframe + pack.mass);
fprintf('Set Ns = %d, Np = %d and cell ''%s'' in the Shepherd block of vtol_config.\n', Ns, Np, cell_name);
