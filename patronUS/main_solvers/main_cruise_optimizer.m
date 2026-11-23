%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% PATRONUS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Performance Analysis, Trajectory Routing
% and Optimized Navigation for Unmanned Systems

%% Main Cruise Segment Optimizer
% By: Alvaro Martinez Blanco (2026)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
%% User version input

versions.AC_version    = "2Wings"; % Aicraft Version: 1Wing, 2Wings
versions.CL_version    = 'analytical'; % CL model to use ("analytical, dummy")
versions.CD_version    = 'analytical'; % CD model to use ("analytical, dummy")
versions.CM_version    = 'analytical'; % CM model to use ("analytical, dummy")


versions.CT_version    = 'analytical'; % CT model to use
versions.CP_version    = 'analytical'; % CP model to use
versions.CH_version    = 'zeros'; % CH model to use

% Battery model: ' ' (ideal energy store, steady cruise only) or 'Shepherd'.
% With 'Shepherd' the steady optimum is followed by the time-dependent
% cruise: maximum range from SOC0 to SOC_f (cruise_config), with the trim
% imposed at every node (2Wings only)
versions.Batt_version = 'Shepherd'; % Battery model to use

% Shepherd only: dynamics and transcription of the time-dependent cruise
% To choose among: "Cruise_2Wings" / "trapezoidal", "hermite_simpson"
versions.dynamics      = "Cruise_2Wings";
versions.transcription = "trapezoidal";

% Versions for constraints
% To choose among: 1Wing, 1Wing_Nacelle, 2Wings
versions.constraints = "2Wings";

% Versions for constraints
% To choose among: "max_range_1W", "max_range_2W"
versions.objective = "max_range_2W";

%% Import Casadi and optimizer

import casadi.*
opti = Opti();

%% Initialize configuration & Generate Models

[params, bounds, models] = cruise_config(versions);

%% Load Registry and Configuration

registry = patronus_registry();

validate_config(versions,registry);

constraint_fcn = registry.constraints(char(versions.constraints));
objective_fcn  = registry.objective(char(versions.objective));




%% Define Decision Variables

% Define Decision variables Dynamically

var_names = fieldnames(bounds.vars);
n_vars = numel(var_names);

% Allocate generic optimization vector

X = opti.variable(n_vars,1);

% Dynamically build structurl mapping and apply bounds and Ini values

X_struct = struct();

vars0 = zeros(n_vars,1);

for i = 1:n_vars
    name = var_names{i};
    X_struct.(name) = X(i);
    v_cfg = bounds.vars.(name);
    opti.subject_to(v_cfg.min <= X(i) <= v_cfg.max);
    vars0(i) = v_cfg.init;
end
%% Apply initial interation for solver
opti.set_initial(X, vars0)

%% Call Objective Function

obj_fun = objective_fcn(X_struct, models, params);
opti.minimize(obj_fun);

%% Call Constraints

[c, ceq,dbg] = constraint_fcn(X_struct, models, params, bounds); % Can have diag as output ( or not)



g_all = [c; ceq];   % keep this handle
n_ineq = numel(c);
n_eq   = numel(ceq);

con_ineq = (c   <= 0);
con_eq   = (ceq == 0);
opti.subject_to(con_ineq);
opti.subject_to(con_eq);

hist = IterHistory();
opti.callback(@(i) hist.record(opti.debug.value(X), opti.debug.value(g_all)));
%% Configure solver (IPOPT)

opts.ipopt.max_iter = 1000;
opts.ipopt.tol = 1e-9; % ConstraintTolerance
opti.solver('ipopt', opts);

%% Run

solve_ok = true; % Status flag

try
    sol_obj    = opti.solve();
    sol_raw    = sol_obj.value(X);
    fval       = sol_obj.value(obj_fun);
    lagmul_ineq = opti.debug.value(opti.dual(con_ineq));
    lagmul_eq   = opti.debug.value(opti.dual(con_eq));
    lam_all     = sol_obj.value(opti.lam_g);

catch e
    solve_ok = false;
    fprintf('Optimization failed or hit limits. Extracting last debug values.\n');
    fprintf('  Reason: %s\n', e.message);
    sol_obj = [];
    sol_raw = opti.debug.value(X);
    fval    = opti.debug.value(obj_fun);
    lam_all = opti.debug.value(opti.lam_g);
    lagmul_ineq = opti.debug.value(opti.dual(con_ineq));
    lagmul_eq   = opti.debug.value(opti.dual(con_eq));
end

if solve_ok

    if versions.AC_version == "1Wing"
    [results, checks, aux] = cruise_postp_1Wing(sol_obj, X_struct, models, params, fval,versions);
    generate_results_page_1Wing(results, models)

    elseif  versions.AC_version == "2Wings"
    if exist('dbg','var')
        [results] = cruise_postp_2Wings_dbg(sol_obj,X_struct,dbg,params,fval);
    else
        [results, checks, aux] = cruise_postp_2Wings(sol_obj, X_struct, models, params, fval,versions);
    end
    generate_results_page_2Wings(results, models,params)
    end
else
    warning('Skipping post-processing: no converged solution.');
end

% IPOPT Diagnosis

stats = opti.stats();
inf_pr = stats.iterations.inf_pr;    % primal infeasibility (constraint violation) per iter
inf_du = stats.iterations.inf_du;    % dual infeasibility (KKT stationarity) per iter
obj_hist = stats.iterations.obj;

viol  = computeConstraintViolation(hist.g_hist, n_ineq, n_eq);
slack = computeBoundSlack(hist.X_hist, bounds);

history_matrix = hist.X_hist;
plotIpoptconvergence(history_matrix, bounds);
plotVariableEvolution(hist.X_hist, bounds);
plotConvergenceMetrics(stats, viol, slack);

%% Time-dependent cruise with the Shepherd battery
% Maximum range for the charge between SOC0 and SOC_f. The steady optimum
% above (ideal-battery objective, min P/V) held constant is the initial
% guess; the optimiser may then change the trim as the voltage drops.

if strcmp(params.batt.model, 'Shepherd')

    if versions.AC_version ~= "2Wings"
        error('main_cruise_optimizer:shepherd', 'The time-dependent cruise is only set up for 2Wings.');
    end

    fns.dynamics    = registry.dynamics(char(versions.dynamics));
    fns.constraints = @(Z, m, p, b) cruise_path_constraints(constraint_fcn, Z, m, p, b);
    fns.objective   = registry.objective('max_range_soc');
    method          = registry.transcription(char(versions.transcription));

    if solve_ok, val = @(e) full(sol_obj.value(e)); else, val = @(e) full(opti.debug.value(e)); end
    trim = struct();
    for i = 1:n_vars
        trim.(var_names{i}) = val(X_struct.(var_names{i}));
    end
    P_bus_steady = val(dbg.P1 + dbg.P2)/params.prop.eff;
    guess = cruise_soc_guess(trim, P_bus_steady, params, bounds);

    [res_soc, info_soc] = ocp_ipopt_stage(params, bounds, models, fns, method, guess, 5);

    fprintf('\nTime-dependent cruise (Shepherd, SOC %.2f -> %.2f): %s, %d iterations, %.1f s\n', ...
            params.batt.SOC0, params.cruise.SOC_f, info_soc.status, info_soc.iterations, info_soc.time);
    fprintf('  Range: %.2f km in %.1f min (steady trim held: %.2f km in %.1f min, %+.2f %%)\n', ...
            res_soc.x(end)/1e3, res_soc.tf/60, guess.x(end)/1e3, guess.tf/60, ...
            100*(res_soc.x(end)/guess.x(end) - 1));
    fprintf('  V: %.2f -> %.2f m/s   P_bus: %.0f -> %.0f W   U: %.1f -> %.1f V   I: %.1f -> %.1f A\n', ...
            res_soc.V(1), res_soc.V(end), res_soc.dbg.P_bus(1), res_soc.dbg.P_bus(end), ...
            res_soc.dbg.U_batt(1), res_soc.dbg.U_batt(end), res_soc.dbg.I_batt(1), res_soc.dbg.I_batt(end));
    fprintf('  Bus energy drawn: %.0f Wh; ideal battery with that energy at the steady optimum: %.2f km\n', ...
            res_soc.E(end)/3600, res_soc.E(end)*trim.V/P_bus_steady/1e3);

    cruise_plot_soc(res_soc, params, guess)
    vtolto_plot_battery(res_soc, params, 'cruise_battery')

end