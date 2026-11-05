function activity = report_constraint_activity(value_fcn, opti, vars, model, ocp, con_ineq, bounds, results)
%REPORT_CONSTRAINT_ACTIVITY Multipliers and slack of the inequality constraints at the solution.
%   Each constraint of the OCP is one NLP constraint per node (and per
%   midpoint for bounds on Hermite-Simpson midpoint controls), each with
%   its own multiplier. For every path constraint row, final-state
%   inequality, actuator rate limit and variable bound (lower and upper
%   side separately) this
%   keeps the per-node arrays and summarises them:
%     t, value, limit: the constrained quantity at each node and its limit.
%                Bounds and actuator rates (at interval centres) in
%                physical units; path rows and final-state inequalities
%                as the normalised constraint c (limit 0).
%     slack:     distance to the limit, normalised (path rows by their
%                limit, bounds by the variable range).
%     lam:       multiplier per node: sensitivity of the objective to the
%                limit at that node. A node's multiplier is roughly the
%                continuous multiplier times the interval length, so it
%                shrinks with N.
%     lam_sum:   sum over nodes, the sensitivity to relaxing the limit along
%                the whole trajectory (independent of N). Rows are sorted by it.
%     n_active:  nodes where lam > slack. IPOPT ends with lam*slack ~ mu
%                (~1e-9) for every inequality, so the larger factor tells
%                whether the constraint binds (lam) or has slack (slack).
%   Multipliers are per unit of the normalised limit; for bounds, per unit
%   of the nondimensional variable (v/scale).
%   Rows with an active node or within 5 % of their limit are printed.
%   plot_constraint_activity draws the per-node arrays.
%   value_fcn: @(e) sol_obj.value(e) or @(e) opti.debug.value(e).

t  = results.t;
N1 = vars.N + 1;

rows = struct([]);

%% Path constraints, final-state inequalities and actuator rates
% con_ineq = [c_path(:); bc_ineq; rate limits]

lam_c    = full(value_fcn(opti.dual(con_ineq)));
lam_path = reshape(lam_c(1:model.n_c*N1), model.n_c, N1);
lam_bc   = lam_c(model.n_c*N1 + (1:ocp.n_bc_ineq));
lam_rate = lam_c(model.n_c*N1 + ocp.n_bc_ineq + 1:end);

for k = 1:model.n_c
    c = results.c_path(k,:);
    rows = [rows, make_row(model.c_names{k}, 'path', t, c, 0, -c, lam_path(k,:))]; %#ok<AGROW>
end
for k = 1:numel(ocp.bc_ineq_names)
    c = full(value_fcn(ocp.bc_ineq(k)));
    rows = [rows, make_row(ocp.bc_ineq_names{k}, 'final', t(end), c, 0, -c, lam_bc(k))]; %#ok<AGROW>
end
for r = ocp.rate_rows
    rate  = full(value_fcn(r.rate));          % Physical rate between control points
    slack = 1 - rate/r.limit;                 % Normalised by the limit
    rows  = [rows, make_row(r.name, 'rate', r.tau*results.tf, rate, r.limit, slack, lam_rate(r.idx)')]; %#ok<AGROW>
end

%% Variable bounds
% Opti splits min <= v <= max into two constraints: the dual is
% [lower bound multipliers, upper bound multipliers], both >= 0

for i = 1:vars.nx
    name  = vars.state_names{i};
    nodes = vars.con.X_nodes{i}; % Nodes fixed by boundary conditions are not bounded
    [lo, up] = bound_dual(value_fcn, opti, vars.con.X{i});
    rows = [rows, bound_rows(name, t(nodes), results.(name)(nodes), lo, up, bounds.states.(name))]; %#ok<AGROW>
end
for i = 1:vars.nu
    name = vars.control_names{i};
    tt   = t;
    vals = results.(name);
    [lo, up] = bound_dual(value_fcn, opti, vars.con.U{i});
    if size(vars.Umn, 2) > 0
        [lo_m, up_m] = bound_dual(value_fcn, opti, vars.con.Um{i});
        [tt, order] = sort([tt, results.t_mid]); % Nodes and midpoints in time order
        vals = [vals, results.mid.(name)];  vals = vals(order);
        lo   = [lo, lo_m];                  lo   = lo(order);
        up   = [up, up_m];                  up   = up(order);
    end
    rows = [rows, bound_rows(name, tt, vals, lo, up, bounds.controls.(name))]; %#ok<AGROW>
end
[lo, up] = bound_dual(value_fcn, opti, vars.con.tf);
rows = [rows, bound_rows('tf', results.tf, results.tf, lo, up, bounds.time.tf)];

%% Report

[~, order] = sort([rows.lam_sum], 'descend');
activity = rows(order);

show = [activity.n_active] > 0 | [activity.slack_min] < 0.05;
fprintf('\nConstraint activity (multiplier = objective decrease per unit relaxation of the normalised limit)\n');
fprintf('  %-22s %10s %11s %11s %8s %11s %8s\n', 'constraint', 'active', 'sum lambda', 'max lambda', 't [s]', 'min slack', 't [s]');
for r = activity(show)
    fprintf('  %-22s %4d/%-5d %11.3e %11.3e %8.3f %11.3e %8.3f\n', r.name, r.n_active, r.n_points, ...
            r.lam_sum, r.lam_max, r.t_lam, r.slack_min, r.t_slack);
end
fprintf('  (%d constraints with no active node and > 5%% slack not shown)\n\n', sum(~show));

end

function [lo, up] = bound_dual(value_fcn, opti, con)
lam = full(value_fcn(opti.dual(con)));
n   = numel(lam)/2;
lo  = reshape(lam(1:n),     1, []);
up  = reshape(lam(n+1:end), 1, []);
end

function rows = bound_rows(name, t, vals, lam_lo, lam_up, cfg)
span = max(cfg.max - cfg.min, eps);
rows = [make_row(sprintf('%s >= %g', name, cfg.min), 'bound', t, vals, cfg.min, (vals - cfg.min)/span, lam_lo), ...
        make_row(sprintf('%s <= %g', name, cfg.max), 'bound', t, vals, cfg.max, (cfg.max - vals)/span, lam_up)];
end

function r = make_row(name, kind, t, value, limit, slack, lam)
lam = abs(lam);
[lam_max,   i] = max(lam);
[slack_min, j] = min(slack);
r = struct('name', name, 'kind', kind, ...
           'n_active', sum(lam > slack), 'n_points', numel(slack), ...
           'lam_sum', sum(lam), 'lam_max', lam_max, 't_lam', t(i), ...
           'slack_min', slack_min, 't_slack', t(j), ...
           't', t, 'value', value, 'limit', limit, 'slack', slack, 'lam', lam);
end
