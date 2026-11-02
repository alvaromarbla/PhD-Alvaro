function obj_fun = vtolto_objective_min_energy(traj, models, params)
% Minimum battery energy consumed during the take-off, in % of E_batt.
% traj holds the node histories by name (see transcribe_ocp).

obj_fun = 100*traj.E(end)/params.E_batt;

end
