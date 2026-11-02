function chk = check_ocp_guess(guess, dynamics_fcn, constraint_fcn, bounds, models, params)
%CHECK_OCP_GUESS Feasibility of an initial guess against the continuous OCP.
%   guess has the layout of unpack_ocp_solution (t, tf, <state>, <control>).
%   Checks, at every point of the guess:
%     bound_viol: state/control box bounds and tf bounds, as a fraction of
%                 each variable's range
%     bc_viol:    boundary conditions x0, xf (and xf_min/xf_max), same scaling
%     path_viol:  max(c) and max|ceq| of constraint_fcn (raw units)
%     dyn_res:    max over states of |f(x,u) - dx/dt| / (x_scale/tf), with
%                 dx/dt the numerical derivative of the guess and x_scale
%                 the state scaling of define_ocp_variables
%   chk.feasible is true when the first three are ~0 and dyn_res is below
%   params.ocp.guess_dyn_tol.

tol = 1e-6;

state_names   = fieldnames(bounds.states);
control_names = fieldnames(bounds.controls);
n_pts = numel(guess.t);

%% Box bounds and boundary conditions

bound_viol = range_violation(guess.tf, bounds.time.tf);
bc_viol    = 0;

for i = 1:numel(state_names)
    cfg = bounds.states.(state_names{i});
    y   = guess.(state_names{i});
    span = max(cfg.max - cfg.min, eps);
    bound_viol = max(bound_viol, range_violation(y, cfg));
    if ~isnan(cfg.x0), bc_viol = max(bc_viol, abs(y(1)   - cfg.x0)/span); end
    if ~isnan(cfg.xf), bc_viol = max(bc_viol, abs(y(end) - cfg.xf)/span); end
    if isfield(cfg, 'xf_max'), bc_viol = max(bc_viol, (y(end) - cfg.xf_max)/span); end
    if isfield(cfg, 'xf_min'), bc_viol = max(bc_viol, (cfg.xf_min - y(end))/span); end
end
for i = 1:numel(control_names)
    bound_viol = max(bound_viol, range_violation(guess.(control_names{i}), bounds.controls.(control_names{i})));
end

%% Path constraints and dynamics, point by point

path_viol = 0;
xdot_mod = zeros(numel(state_names), n_pts);

for k = 1:n_pts
    Z = struct();
    for i = 1:numel(state_names)
        Z.(state_names{i}) = guess.(state_names{i})(k);
    end
    for i = 1:numel(control_names)
        Z.(control_names{i}) = guess.(control_names{i})(k);
    end

    [c, ceq] = constraint_fcn(Z, models, params, bounds);
    path_viol = max([path_viol; full(c(:)); abs(full(ceq(:)))]);

    xdot = dynamics_fcn(Z, models, params);
    for i = 1:numel(state_names)
        xdot_mod(i,k) = full(xdot.(state_names{i}));
    end
end

dyn_res = 0;
for i = 1:numel(state_names)
    cfg = bounds.states.(state_names{i});
    xdot_num = derivative_2nd_order(guess.(state_names{i}), guess.t);
    scale    = max([abs(cfg.min), abs(cfg.max), 1])/guess.tf; % Same state scaling as the defects
    dyn_res  = max(dyn_res, max(abs(xdot_mod(i,:) - xdot_num))/scale);
end

%% Output

chk.bound_viol = bound_viol;
chk.bc_viol    = bc_viol;
chk.path_viol  = path_viol;
chk.dyn_res    = dyn_res;
chk.feasible   = bound_viol <= tol && bc_viol <= tol && path_viol <= tol ...
                 && dyn_res <= params.ocp.guess_dyn_tol;

end

function dy = derivative_2nd_order(y, t)
% Central differences inside, 2nd-order one-sided at the ends (uniform t).
% gradient() is 1st order at the ends, which shows up as a false residual
% where the guess has a large jerk (start/end of the splines)
dt = t(2) - t(1);
dy = gradient(y, dt);
dy(1)   = (-3*y(1)   + 4*y(2)     - y(3))/(2*dt);
dy(end) = ( 3*y(end) - 4*y(end-1) + y(end-2))/(2*dt);
end

function v = range_violation(y, cfg)
% Largest excursion outside [min, max], as a fraction of the range
span = max(cfg.max - cfg.min, eps);
v = max([0, (cfg.min - min(y))/span, (max(y) - cfg.max)/span]);
end
