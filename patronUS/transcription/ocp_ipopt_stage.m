function [results, info] = ocp_ipopt_stage(params, bounds, models, fns, method, guess, opts)
%OCP_IPOPT_STAGE One transcription stage of a trajectory segment solved with IPOPT.
%   Same building blocks as main_vtolto_optimizer, as a function, for
%   segments and studies that solve an OCP in one call (time-dependent
%   cruise, main_battery_sizing).
%   fns: struct with dynamics, constraints and objective function handles
%   method: entry of registry.transcription; guess: trajectory struct or []
%   opts: print_level (IPOPT print level, default 0 = silent),
%         activity (default false: true reports the multipliers and slack
%         of every inequality and bound, report_constraint_activity, in
%         info.activity for plot_constraint_activity)
%   results: unpack_ocp_solution output; info: ok, status, objective,
%   iterations, time (and activity)

import casadi.*

if nargin < 7, opts = struct(); end
if ~isfield(opts, 'print_level'), opts.print_level = 0;     end
if ~isfield(opts, 'activity'),    opts.activity    = false; end

opti = Opti();
vars = define_ocp_variables(opti, bounds, params.ocp.N, method);
if ~isempty(guess), set_ocp_initial_guess(opti, vars, guess); end

model = build_continuous_model(fns.dynamics, fns.constraints, vars, models, params, bounds);
ocp   = transcribe_ocp(vars, model, method, bounds);
obj   = fns.objective(ocp.traj, models, params) + control_rate_penalty(vars, params.ocp.reg_weight);
opti.minimize(obj);
con_ineq = (ocp.c <= 0);
opti.subject_to(con_ineq);
opti.subject_to(ocp.ceq == 0);
opti.solver('ipopt', struct('print_time', false), struct('max_iter', 3000, 'tol', 1e-8, 'print_level', opts.print_level));

t0 = tic;
try
    sol = opti.solve();
    info.ok = true;
    value_fcn = @(e) sol.value(e);
catch
    info.ok = false;
    value_fcn = @(e) opti.debug.value(e);
end
info.time       = toc(t0);
info.iterations = opti.stats().iter_count;
info.status     = opti.stats().return_status;
info.objective  = full(value_fcn(obj));

results = unpack_ocp_solution(value_fcn, vars, model);

if opts.activity
    info.activity = report_constraint_activity(value_fcn, opti, vars, model, ocp, con_ineq, bounds, results);
end

end
