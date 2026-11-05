function ocp = transcribe_ocp(vars, model, method, bounds)
%TRANSCRIBE_OCP Transcribe the continuous OCP into constraints on the discrete variables.
%   On the uniform grid t_k = tau_k*tf, h = tf/N, it builds:
%     - collocation defects with method.defects (trapezoidal, Hermite-Simpson, ...)
%     - path constraints of the segment at every node
%     - boundary conditions x(0) = x0, x(tf) = xf from bounds.states (NaN = free)
%       and, if present, xf_min <= x(tf) <= xf_max
%     - terminal derivative conditions dx/dt(tf) = dxf, if present: with
%       dxf = 0 on velocities and rates the final state is an equilibrium
%       (trimmed), not just a point where the velocity crosses zero
%   Defects and boundary conditions are divided by vars.x_scale so all
%   states weigh alike in IPOPT. Output follows the cruise convention:
%     ocp.c <= 0,  ocp.ceq == 0  with ceq = [defects; path ceq; boundary]
%   ocp.traj holds the symbolic node histories by name (plus t, tf) and is
%   what the objective functions receive.

N = vars.N;
h = vars.tf/N;

%% Continuous model at the nodes

F_nodes = model.F.map(N+1);
[Fk, c_path, ceq_path] = F_nodes(vars.X, vars.U);

%% Collocation defects

defects = method.defects(model.F, vars, Fk, h);
defects = defects ./ repmat(vars.x_scale, 1, N);

%% Boundary conditions

bc      = {};
bc_ineq = {};
bc_ineq_names = {};
for i = 1:vars.nx
    name = vars.state_names{i};
    cfg  = bounds.states.(name);
    if ~isnan(cfg.x0)
        bc{end+1, 1} = (vars.X(i,1)   - cfg.x0)/vars.x_scale(i); %#ok<AGROW>
    end
    if ~isnan(cfg.xf)
        bc{end+1, 1} = (vars.X(i,end) - cfg.xf)/vars.x_scale(i); %#ok<AGROW>
    end
    if isfield(cfg, 'xf_max')
        bc_ineq{end+1, 1} = (vars.X(i,end) - cfg.xf_max)/vars.x_scale(i); %#ok<AGROW>
        bc_ineq_names{end+1} = sprintf('%s(tf) <= %g', name, cfg.xf_max); %#ok<AGROW>
    end
    if isfield(cfg, 'xf_min')
        bc_ineq{end+1, 1} = (cfg.xf_min - vars.X(i,end))/vars.x_scale(i); %#ok<AGROW>
        bc_ineq_names{end+1} = sprintf('%s(tf) >= %g', name, cfg.xf_min); %#ok<AGROW>
    end
    if isfield(cfg, 'dxf') && ~isnan(cfg.dxf)
        bc{end+1, 1} = (Fk(i,end) - cfg.dxf)/vars.x_scale(i); %#ok<AGROW>
    end
end
bc      = vertcat(bc{:});
bc_ineq = vertcat(bc_ineq{:});

%% NLP constraints

ocp.c   = [c_path(:); bc_ineq];
ocp.ceq = [defects(:); ceq_path(:); bc];

ocp.n_defects = vars.nx*N;
ocp.n_path_c  = model.n_c*(N+1);
ocp.n_path_eq = model.n_ceq*(N+1);
ocp.n_bc      = size(bc, 1);
ocp.n_bc_ineq = size(bc_ineq, 1);
ocp.bc_ineq       = bc_ineq;       % For report_constraint_activity
ocp.bc_ineq_names = bc_ineq_names;
ocp.n_ineq    = ocp.n_path_c + ocp.n_bc_ineq;
ocp.n_eq      = ocp.n_defects + ocp.n_path_eq + ocp.n_bc;

%% Symbolic trajectory by name (for objectives)

ocp.traj.tf = vars.tf;
ocp.traj.t  = vars.tau*vars.tf;
for i = 1:vars.nx
    ocp.traj.(vars.state_names{i}) = vars.X(i,:);
end
for i = 1:vars.nu
    ocp.traj.(vars.control_names{i}) = vars.U(i,:);
end

end
