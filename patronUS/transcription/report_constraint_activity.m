function activity = report_constraint_activity(value_fcn, opti, vars, model, ocp, con_ineq, bounds, results)
%REPORT_CONSTRAINT_ACTIVITY Multipliers and slack of the inequality constraints at the solution.
%   For every path constraint row, final-state inequality and variable
%   bound (lower and upper side separately):
%     n_active:  points where the constraint is active (slack <= 1e-6)
%     lam_max:   largest |multiplier| and the time t_lam where it occurs.
%                The multiplier is the sensitivity of the objective to the
%                limit: relaxing it by one (normalised) unit lowers the
%                objective by about lam. Zero means the constraint does not
%                shape the solution, even if it is touched. For bounds the
%                unit is the nondimensional variable (v/scale).
%     slack_min: smallest distance to the limit and its time t_slack.
%                Path rows are normalised by their limit (vtolto path
%                constraints), bounds by the variable range.
%   Rows with a multiplier or within 5 % of their limit are printed, sorted
%   by lam_max. activity holds all rows.
%   value_fcn: @(e) sol_obj.value(e) or @(e) opti.debug.value(e).

tol = 1e-6;
t     = results.t;
N1    = vars.N + 1;

rows = struct('name', {}, 'n_active', {}, 'n_points', {}, 'lam_max', {}, ...
              't_lam', {}, 'slack_min', {}, 't_slack', {});

%% Path constraints and final-state inequalities (con_ineq = [c_path(:); bc_ineq])

lam_c    = full(value_fcn(opti.dual(con_ineq)));
lam_path = reshape(lam_c(1:model.n_c*N1), model.n_c, N1);
lam_bc   = lam_c(model.n_c*N1+1:end);

for k = 1:model.n_c
    rows(end+1) = make_row(model.c_names{k}, -results.c_path(k,:), lam_path(k,:), t, tol); %#ok<AGROW>
end
for k = 1:numel(ocp.bc_ineq_names)
    slack = -full(value_fcn(ocp.bc_ineq(k)));
    rows(end+1) = make_row(ocp.bc_ineq_names{k}, slack, lam_bc(k), t(end), tol); %#ok<AGROW>
end

%% Variable bounds
% Opti splits min <= v <= max into two constraints: the dual is
% [lower bound multipliers, upper bound multipliers], both >= 0

for i = 1:vars.nx
    name = vars.state_names{i};
    [lo, up] = bound_dual(value_fcn, opti, vars.con.X{i});
    rows = [rows, bound_rows(name, results.(name), lo, up, t, bounds.states.(name), tol)]; %#ok<AGROW>
end
for i = 1:vars.nu
    name = vars.control_names{i};
    vals = results.(name);
    [lo, up] = bound_dual(value_fcn, opti, vars.con.U{i});
    tt   = t;
    if size(vars.Umn, 2) > 0
        [lo_m, up_m] = bound_dual(value_fcn, opti, vars.con.Um{i});
        vals = [vals, results.mid.(name)];
        lo   = [lo, lo_m];
        up   = [up, up_m];
        tt   = [tt, results.t_mid];
    end
    rows = [rows, bound_rows(name, vals, lo, up, tt, bounds.controls.(name), tol)]; %#ok<AGROW>
end
[lo, up] = bound_dual(value_fcn, opti, vars.con.tf);
rows = [rows, bound_rows('tf', results.tf, lo, up, results.tf, bounds.time.tf, tol)];

%% Report

[~, order] = sort([rows.lam_max], 'descend');
activity = rows(order);

show = [activity.lam_max] > 1e-8 | [activity.slack_min] < 0.05;
fprintf('\nConstraint activity (multiplier = objective decrease per unit relaxation of the normalised limit)\n');
fprintf('  %-22s %10s %11s %8s %11s %8s\n', 'constraint', 'active', 'max|lambda|', 't [s]', 'min slack', 't [s]');
for r = activity(show)
    fprintf('  %-22s %4d/%-5d %11.3e %8.3f %11.3e %8.3f\n', r.name, r.n_active, r.n_points, ...
            r.lam_max, r.t_lam, r.slack_min, r.t_slack);
end
fprintf('  (%d constraints with no multiplier and > 5%% slack not shown)\n\n', sum(~show));

end

function [lo, up] = bound_dual(value_fcn, opti, con)
lam = full(value_fcn(opti.dual(con)));
n   = numel(lam)/2;
lo  = reshape(lam(1:n),     1, []);
up  = reshape(lam(n+1:end), 1, []);
end

function rows = bound_rows(name, vals, lam_lo, lam_up, t, cfg, tol)
span = max(cfg.max - cfg.min, eps);
rows = [make_row(sprintf('%s >= %g', name, cfg.min), (vals - cfg.min)/span, lam_lo, t, tol), ...
        make_row(sprintf('%s <= %g', name, cfg.max), (cfg.max - vals)/span, lam_up, t, tol)];
end

function r = make_row(name, slack, lam, t, tol)
[lam_max,   i] = max(abs(lam));
[slack_min, j] = min(slack);
r = struct('name', name, 'n_active', sum(slack <= tol), 'n_points', numel(slack), ...
           'lam_max', lam_max, 't_lam', t(i), 'slack_min', slack_min, 't_slack', t(j));
end
