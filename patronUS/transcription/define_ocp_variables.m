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

import casadi.*

vars.state_names   = fieldnames(bounds.states);
vars.control_names = fieldnames(bounds.controls);
vars.nx = numel(vars.state_names);
vars.nu = numel(vars.control_names);
vars.N  = N;

vars.tau     = linspace(0, 1, N+1);
vars.tau_mid = 0.5*(vars.tau(1:end-1) + vars.tau(2:end));

%% Allocate variables

vars.X  = opti.variable(vars.nx, N+1);
vars.U  = opti.variable(vars.nu, N+1);
if method.midpoints
    vars.Um = opti.variable(vars.nu, N);
else
    vars.Um = MX(vars.nu, 0);
end
vars.tf = opti.variable();

%% States: bounds, initial guess and scaling

vars.x_scale = ones(vars.nx, 1); % Used to scale defects and boundary conditions

for i = 1:vars.nx
    cfg = bounds.states.(vars.state_names{i});
    opti.subject_to(cfg.min <= vars.X(i,:) <= cfg.max);
    opti.set_initial(vars.X(i,:), guess_profile(cfg.init, vars.tau));
    vars.x_scale(i) = max([abs(cfg.min), abs(cfg.max), 1]);
end

%% Controls: bounds and initial guess

for i = 1:vars.nu
    cfg = bounds.controls.(vars.control_names{i});
    opti.subject_to(cfg.min <= vars.U(i,:) <= cfg.max);
    opti.set_initial(vars.U(i,:), guess_profile(cfg.init, vars.tau));
    if method.midpoints
        opti.subject_to(cfg.min <= vars.Um(i,:) <= cfg.max);
        opti.set_initial(vars.Um(i,:), guess_profile(cfg.init, vars.tau_mid));
    end
end

%% Final time

cfg = bounds.time.tf;
opti.subject_to(cfg.min <= vars.tf <= cfg.max);
opti.set_initial(vars.tf, cfg.init);

%% Flat decision vector (iteration history / bound slack)
% Layout: [X(:); U(:); Um(:); tf], column-major (node by node)

vars.w = [vars.X(:); vars.U(:); vars.Um(:); vars.tf];

end

function y = guess_profile(init, tau)
% Constant guess for a scalar, linear in tau for [start end]
if isscalar(init)
    y = init*ones(size(tau));
else
    y = init(1) + (init(end) - init(1))*tau;
end
end
