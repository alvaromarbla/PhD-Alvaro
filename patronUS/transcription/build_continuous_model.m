function model = build_continuous_model(dynamics_fcn, constraint_fcn, vars, models, params, bounds)
%BUILD_CONTINUOUS_MODEL Wrap the continuous dynamics and path constraints in one CasADi Function.
%   model.F(x, u) -> [xdot, c, ceq, dbg] for a single node. It is built once
%   from the name-based segment functions and then mapped over the grid by
%   the collocation schemes, which keeps the NLP graph small.
%
%   dynamics_fcn:   [xdot, dbg] = f(Z, models, params),  xdot.<state_name>
%   constraint_fcn: [c, ceq, dbg] = f(Z, models, params, bounds), c <= 0, ceq == 0
%   Z holds one scalar per state and control name (like X_struct in cruise).
%   Scalar fields of both dbg structs are stacked into the 'dbg' output so
%   they can be evaluated at the solution (see unpack_ocp_solution).

import casadi.*

x = MX.sym('x', vars.nx);
u = MX.sym('u', vars.nu);

Z = struct();
for i = 1:vars.nx
    Z.(vars.state_names{i}) = x(i);
end
for i = 1:vars.nu
    Z.(vars.control_names{i}) = u(i);
end

[xdot_s, dbg_dyn]   = dynamics_fcn(Z, models, params);
[c, ceq, dbg_con]   = constraint_fcn(Z, models, params, bounds);

%% State derivatives, ordered as the rows of X

missing = setdiff(vars.state_names, fieldnames(xdot_s));
if ~isempty(missing)
    error('build_continuous_model:missingDerivative', ...
        'Dynamics do not return a derivative for state(s): %s', strjoin(missing, ', '));
end
xdot = cellfun(@(name) xdot_s.(name), vars.state_names, 'UniformOutput', false);
xdot = vertcat(xdot{:});

%% Path constraints

if isempty(c),   c   = MX(0,1); end
if isempty(ceq), ceq = MX(0,1); end

%% Debug outputs (scalar fields only, path constraint fields override dynamics ones)

dbg = dbg_dyn;
con_names = fieldnames(dbg_con);
for k = 1:numel(con_names)
    dbg.(con_names{k}) = dbg_con.(con_names{k});
end
dbg_names = fieldnames(dbg);
dbg_names = dbg_names(cellfun(@(name) isequal(size(dbg.(name)), [1 1]), dbg_names));
dbg_vec   = cellfun(@(name) dbg.(name), dbg_names, 'UniformOutput', false);
dbg_vec   = vertcat(MX(0,1), dbg_vec{:});

%% CasADi Function

model.F = Function('F_continuous', {x, u}, {xdot, c, ceq, dbg_vec}, ...
                   {'x', 'u'}, {'xdot', 'c', 'ceq', 'dbg'});
model.n_c       = size(c, 1);
model.n_ceq     = size(ceq, 1);
model.dbg_names = dbg_names;

end
