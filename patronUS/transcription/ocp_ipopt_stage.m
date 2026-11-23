function [results, info] = ocp_ipopt_stage(params, bounds, models, fns, method, guess, print_level)
%OCP_IPOPT_STAGE One transcription stage of a trajectory segment solved with IPOPT.
%   Same building blocks as main_vtolto_optimizer, as a function, for
%   segments and studies that solve an OCP in one call (time-dependent
%   cruise, main_battery_sizing).
%   fns: struct with dynamics, constraints and objective function handles
%   method: entry of registry.transcription; guess: trajectory struct or []
%   print_level: IPOPT print level (default 0, silent)
%   results: unpack_ocp_solution output; info: ok, objective, iterations, time

import casadi.*

if nargin < 7, print_level = 0; end

opti = Opti();
vars = define_ocp_variables(opti, bounds, params.ocp.N, method);
if ~isempty(guess), set_ocp_initial_guess(opti, vars, guess); end

model = build_continuous_model(fns.dynamics, fns.constraints, vars, models, params, bounds);
ocp   = transcribe_ocp(vars, model, method, bounds);
obj   = fns.objective(ocp.traj, models, params) + control_rate_penalty(vars, params.ocp.reg_weight);
opti.minimize(obj);
opti.subject_to(ocp.c <= 0);
opti.subject_to(ocp.ceq == 0);
opti.solver('ipopt', struct('print_time', false), struct('max_iter', 3000, 'tol', 1e-8, 'print_level', print_level));

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

end
