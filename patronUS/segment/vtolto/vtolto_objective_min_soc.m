function obj_fun = vtolto_objective_min_soc(traj, models, params)
% Minimum state-of-charge drop during the take-off, in % (Shepherd battery).
% Unlike the energy objective, it sees the voltage sag: at high power the
% battery voltage drops and the same energy costs more charge.

if ~isfield(traj, 'SOC')
    error('vtolto_objective_min_soc:battery', 'min_soc_vtol needs versions.Batt_version = ''Shepherd''.');
end
obj_fun = 100*(traj.SOC(1) - traj.SOC(end));

end
