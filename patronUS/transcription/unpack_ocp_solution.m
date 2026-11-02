function res = unpack_ocp_solution(value_fcn, vars, model)
%UNPACK_OCP_SOLUTION Numerical trajectory at the nodes, by name.
%   value_fcn evaluates a symbolic expression at the solution, e.g.
%   @(e) sol_obj.value(e) or @(e) opti.debug.value(e) on a failed solve.
%   res.t, res.tf, res.<state>, res.<control> are 1 x (N+1) rows;
%   res.xdot.<state> are the state derivatives and res.dbg.<name> the debug
%   expressions of the segment functions, evaluated at every node.

X  = value_fcn(vars.X);
U  = value_fcn(vars.U);
tf = value_fcn(vars.tf);

res.tf = tf;
res.t  = vars.tau*tf;

for i = 1:vars.nx
    res.(vars.state_names{i}) = X(i,:);
end
for i = 1:vars.nu
    res.(vars.control_names{i}) = U(i,:);
end
if size(vars.Um, 2) > 0
    Um = value_fcn(vars.Um);
    res.t_mid = vars.tau_mid*tf;
    for i = 1:vars.nu
        res.mid.(vars.control_names{i}) = Um(i,:);
    end
end

%% Re-evaluate the continuous model at the solution nodes

F_nodes = model.F.map(vars.N+1);
[xdot, c, ceq, dbg] = F_nodes(X, U);

xdot = full(xdot);
dbg  = full(dbg);
for i = 1:vars.nx
    res.xdot.(vars.state_names{i}) = xdot(i,:);
end
for k = 1:numel(model.dbg_names)
    res.dbg.(model.dbg_names{k}) = dbg(k,:);
end
res.c_path   = full(c);
res.ceq_path = full(ceq);

end
