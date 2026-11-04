function vars = define_ocp_variables(opti, bounds, N, method)
%DEFINE_OCP_VARIABLES Discrete decision variables of a direct-collocation OCP.
%   States live at the N+1 nodes of a uniform grid in tau = t/tf in [0, 1],
%   controls at the nodes (and at the N interval midpoints when
%   method.midpoints is true, e.g. Hermite-Simpson), plus the final time tf.
%
%   Each field of bounds.states / bounds.controls defines one row of X / U,
%   in fieldnames order, so adding/commenting a field in the config changes
%   the problem (same idea as bounds.vars in cruise).
%     bounds.states.<name>   = struct('min','max','x0','xf','init'[,'xf_min','xf_max'])
%     bounds.controls.<name> = struct('min','max','init')
%     bounds.time.tf         = struct('min','max','init')
%   'init' is the initial guess: a scalar (constant) or [start end] (linear
%   in time). Box bounds and initial guess are applied here; the boundary
%   conditions x0/xf are imposed in transcribe_ocp.
%
%   Nondimensional variables: IPOPT works on Xn, Un, Umn = X, U, Um divided
%   by x_scale / u_scale = max(|min|, |max|, 1) of each row, so all of them
%   are O(1). vars.X, vars.U, vars.Um are the physical expressions used by
%   the rest of the code; initial guesses are set through set_ocp_initial_guess.
%   The bound constraints are kept in vars.con so their multipliers can be
%   read after the solve (report_constraint_activity).

import casadi.*

vars.state_names   = fieldnames(bounds.states);
vars.control_names = fieldnames(bounds.controls);
vars.nx = numel(vars.state_names);
vars.nu = numel(vars.control_names);
vars.N  = N;

vars.tau     = linspace(0, 1, N+1);
vars.tau_mid = 0.5*(vars.tau(1:end-1) + vars.tau(2:end));

%% Scaling

vars.x_scale = cellfun(@(name) row_scale(bounds.states.(name)),   vars.state_names);
vars.u_scale = cellfun(@(name) row_scale(bounds.controls.(name)), vars.control_names);

%% Allocate nondimensional variables and physical expressions

vars.Xn = opti.variable(vars.nx, N+1);
vars.Un = opti.variable(vars.nu, N+1);
if method.midpoints
    vars.Umn = opti.variable(vars.nu, N);
else
    vars.Umn = MX(vars.nu, 0);
end
vars.tf = opti.variable();

vars.X  = diag(vars.x_scale)*vars.Xn;
vars.U  = diag(vars.u_scale)*vars.Un;
vars.Um = diag(vars.u_scale)*vars.Umn;

%% States: bounds and initial guess
% No box bound at the nodes fixed by a boundary condition (x0 at the first
% node, xf at the last): there the bound duplicates the equality, which
% makes the constraint Jacobian degenerate (vars.con.X_nodes keeps the
% bounded nodes)

for i = 1:vars.nx
    cfg = bounds.states.(vars.state_names{i});
    s   = vars.x_scale(i);
    nodes = 1:N+1;
    if ~isnan(cfg.x0), nodes = nodes(nodes ~= 1);   end
    if ~isnan(cfg.xf), nodes = nodes(nodes ~= N+1); end
    vars.con.X_nodes{i} = nodes;
    vars.con.X{i} = (cfg.min/s <= vars.Xn(i,nodes) <= cfg.max/s);
    opti.subject_to(vars.con.X{i});
    opti.set_initial(vars.Xn(i,:), guess_profile(cfg.init, vars.tau)/s);
end

%% Controls: bounds and initial guess

vars.con.Um = cell(1, vars.nu);
for i = 1:vars.nu
    cfg = bounds.controls.(vars.control_names{i});
    s   = vars.u_scale(i);
    vars.con.U{i} = (cfg.min/s <= vars.Un(i,:) <= cfg.max/s);
    opti.subject_to(vars.con.U{i});
    opti.set_initial(vars.Un(i,:), guess_profile(cfg.init, vars.tau)/s);
    if method.midpoints
        vars.con.Um{i} = (cfg.min/s <= vars.Umn(i,:) <= cfg.max/s);
        opti.subject_to(vars.con.Um{i});
        opti.set_initial(vars.Umn(i,:), guess_profile(cfg.init, vars.tau_mid)/s);
    end
end

%% Final time

cfg = bounds.time.tf;
vars.con.tf = (cfg.min <= vars.tf <= cfg.max);
opti.subject_to(vars.con.tf);
opti.set_initial(vars.tf, cfg.init);

%% Flat decision vector in physical units (iteration history / bound slack)
% Layout: [X(:); U(:); Um(:); tf], column-major (node by node)

vars.w = [vars.X(:); vars.U(:); vars.Um(:); vars.tf];

end

function s = row_scale(cfg)
s = max([abs(cfg.min), abs(cfg.max), 1]);
end

function y = guess_profile(init, tau)
% Constant guess for a scalar, linear in tau for [start end]
if isscalar(init)
    y = init*ones(size(tau));
else
    y = init(1) + (init(end) - init(1))*tau;
end
end
