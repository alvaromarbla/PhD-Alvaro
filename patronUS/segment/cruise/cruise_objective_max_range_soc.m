function obj = cruise_objective_max_range_soc(traj, models, params) %#ok<INUSD>
% Maximum range of the time-dependent cruise for a given battery charge
% (SOC from SOC0 to SOC_f, fixed by the boundary conditions): minimise
% -x(tf), non-dimensionalised by params.ocp.x_ref.

obj = -traj.x(end)/params.ocp.x_ref;

end
