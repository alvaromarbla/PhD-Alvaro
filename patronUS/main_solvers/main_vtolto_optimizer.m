%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% PATRONUS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Performance Analysis, Trajectory Routing
% and Optimized Navigation for Unmanned Systems

%% Main VTOL Take-Off Segment Optimizer
% By: Alvaro Martinez Blanco (2026)
% Direct collocation: the continuous OCP is transcribed on N intervals
% (transcription/) and the resulting NLP is solved with IPOPT.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
%% User version input

versions.AC_version    = "2Wings"; % Aicraft Version: 2Wings
versions.CL_version    = 'analytical'; % CL model to use ("analytical, dummy")
versions.CD_version    = 'analytical'; % CD model to use ("analytical, dummy")
versions.CM_version    = 'analytical'; % CM model to use ("analytical, dummy")


versions.CT_version    = 'analytical'; % CT model to use
versions.CP_version    = 'analytical'; % CP model to use
versions.CH_version    = 'zeros'; % CH model to use

versions.Batt_version = ' '; % Battery model to use

% Versions for dynamics
% To choose among: "VTOL_2Wings"
versions.dynamics = "VTOL_2Wings";

% Versions for (path) constraints
% To choose among: "VTOL_2Wings"
versions.constraints = "VTOL_2Wings";

% Versions for objective
% To choose among: "min_energy_vtol", "min_time_vtol"
versions.objective = "min_energy_vtol";

% Versions for initial guess (of the first stage)
% To choose among: "config" (guess in vtol_config), "flat_2Wings" (differential flatness)
versions.initial_guess = "flat_2Wings";

% Versions for transcription. Stages solved in sequence, each one
% warm-started with the solution of the previous one.
% To choose among: "trapezoidal", "hermite_simpson"
versions.transcription = ["trapezoidal", "hermite_simpson"];

%% Import Casadi

import casadi.*

%% Initialize configuration & Generate Models

[params, bounds, models] = vtol_config(versions);

%% Load Registry and Configuration

registry = patronus_registry();

validate_config(versions,registry);

constraint_fcn = registry.constraints(char(versions.constraints));
objective_fcn  = registry.objective(char(versions.objective));
dynamics_fcn   = registry.dynamics(char(versions.dynamics));
guess_fcn      = registry.initial_guess(char(versions.initial_guess));

%% Initial guess of the first stage ([] -> guess in vtol_config)

% Built for each tf in params.ocp.tf_sweep; the feasible one (bounds, path
% constraints, dynamics) with the lowest objective is kept
[guess, guess_sweep] = select_initial_guess(guess_fcn, params.ocp.tf_sweep, objective_fcn, ...
                                            dynamics_fcn, constraint_fcn, bounds, models, params);
warm_start = guess;

%% Transcription stages

for stage = 1:numel(versions.transcription)

    method = registry.transcription(char(versions.transcription(stage)));
    fprintf('\n=== Stage %d/%d: %s ===\n', stage, numel(versions.transcription), versions.transcription(stage));

    opti = Opti();

    %% Define Decision Variables

    % States at the N+1 nodes, controls at the nodes (and midpoints if the
    % method needs them) and final time. Each field of bounds.states /
    % bounds.controls defines one row; bounds and config guess applied inside.

    vars = define_ocp_variables(opti, bounds, params.ocp.N, method);

    if ~isempty(warm_start)
        set_ocp_initial_guess(opti, vars, warm_start);
    end

    %% Continuous model: dynamics + path constraints as one CasADi Function

    model = build_continuous_model(dynamics_fcn, constraint_fcn, vars, models, params, bounds);

    %% Transcription: continuous -> discrete (defects, path constraints, boundary conditions)

    ocp = transcribe_ocp(vars, model, method, bounds);

    %% Call Objective Function

    obj_fun = objective_fcn(ocp.traj, models, params);
    opti.minimize(obj_fun);

    %% Call Constraints

    c   = ocp.c;
    ceq = ocp.ceq;

    g_all = [c; ceq];   % keep this handle
    n_ineq = ocp.n_ineq;
    n_eq   = ocp.n_eq;

    con_ineq = (c   <= 0);
    con_eq   = (ceq == 0);
    opti.subject_to(con_ineq);
    opti.subject_to(con_eq);

    hist = IterHistory();
    opti.callback(@(i) hist.record(opti.debug.value(vars.w), opti.debug.value(g_all)));

    %% Configure solver (IPOPT)

    opts.ipopt.max_iter = 3000;
    opts.ipopt.tol = 1e-8; % ConstraintTolerance
    opti.solver('ipopt', opts);

    %% Run

    solve_ok = true; % Status flag

    try
        sol_obj    = opti.solve();
        fval       = sol_obj.value(obj_fun);
        lagmul_ineq = opti.debug.value(opti.dual(con_ineq));
        lagmul_eq   = opti.debug.value(opti.dual(con_eq));
        lam_all     = sol_obj.value(opti.lam_g);
        value_fcn   = @(e) sol_obj.value(e);

    catch e
        solve_ok = false;
        fprintf('Optimization failed or hit limits. Extracting last debug values.\n');
        fprintf('  Reason: %s\n', e.message);
        sol_obj = [];
        fval    = opti.debug.value(obj_fun);
        lam_all = opti.debug.value(opti.lam_g);
        lagmul_ineq = opti.debug.value(opti.dual(con_ineq));
        lagmul_eq   = opti.debug.value(opti.dual(con_eq));
        value_fcn   = @(e) opti.debug.value(e);
    end

    % Trajectory at the nodes, by name (also on a failed solve, for inspection)
    results = unpack_ocp_solution(value_fcn, vars, model);

    if ~solve_ok
        warning('Stage %d (%s) did not converge: later stages skipped.', stage, versions.transcription(stage));
        break
    end
    warm_start = results;

end

%% Post-processing (last stage)

if solve_ok
    vtolto_plot_results(results, params)
else
    warning('Skipping post-processing: no converged solution.');
end

% IPOPT Diagnosis (last stage)

stats = opti.stats();
inf_pr = stats.iterations.inf_pr;    % primal infeasibility (constraint violation) per iter
inf_du = stats.iterations.inf_du;    % dual infeasibility (KKT stationarity) per iter
obj_hist = stats.iterations.obj;

viol  = computeConstraintViolation(hist.g_hist, n_ineq, n_eq);
slack = ocp_bound_slack(hist.X_hist, vars, bounds);

plotConvergenceMetrics(stats, viol, slack);
