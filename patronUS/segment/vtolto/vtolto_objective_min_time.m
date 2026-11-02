function obj_fun = vtolto_objective_min_time(traj, models, params)
% Minimum take-off time, non-dimensionalised with params.ocp.tf_ref.

obj_fun = traj.tf/params.ocp.tf_ref;

end
