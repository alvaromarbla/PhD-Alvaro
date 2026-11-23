function [c, ceq, dbg] = cruise_path_constraints(steady_fcn, Z, models, params, bounds)
% Path constraints of the time-dependent cruise, imposed at every node:
% the steady cruise constraints (trim equalities and propulsive limits of
% steady_fcn, e.g. cruise_constraints_3GDL, unchanged) plus, with the
% Shepherd battery, the same battery limits as the VTOL take-off:
%     disc >= 0      the pack can deliver P_bus (normalised by U_oc^2)
%     U >= V_cut     terminal voltage above the cut-off
% Z holds the trim variables (controls) and the states x, E, SOC by name.
% Used as @(Z, m, p, b) cruise_path_constraints(steady_fcn, Z, m, p, b).

[c, ceq, dbg] = steady_fcn(Z, models, params, bounds);

if strcmp(params.batt.model, 'Shepherd')
    P_bus = (dbg.P1 + dbg.P2)/params.prop.eff;
    batt  = battery_shepherd(P_bus, Z.SOC, params.batt);
    c = [c;
         -batt.disc/batt.U_oc^2;
         1 - batt.U/params.batt.V_cut];
    n_steady = size(c, 1) - 2;
    dbg.c_names = [arrayfun(@(k) sprintf('trim c%d', k), 1:n_steady, 'UniformOutput', false), ...
                   {'battery power', 'U_batt >= V_cut'}];
end

end
