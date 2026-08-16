%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PATRONUS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Performance Analysis, Trajectory Routing
% and Optimized Navigation for Unmanned Systems

%% Main Cruise Segment Optimizer
% By: Alvaro Martinez Blanco (2026)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all 
%% User version input

versions.AC_version   = "2Wings"; % Aicraft Version
versions.CL_version   = 'analytical'; % CL model to use
versions.CD_version   = 'analytical'; % CD model to use
versions.CM_version   = 'analytical'; % CM model to use


versions.CT_version   = 'analytical'; % CT model to use
versions.CP_version   = 'analytical'; % CP model to use
versions.CH_version   = 'zeros'; % CH model to use

versions.Batt_version = ' '; % Battery model to use

% Versions for constraints
% To choose among: 1Wing, 1Wing_Nacelle, 2Wings
versions.constraints = "2Wings";

% Versions for constraints
% To choose among: "max_range_1W", "max_range_2W"
versions.objective = "max_range_2W";

%% Import Casadi and optimizer

import casadi.*
opti = Opti();

%% Load Registry and Configuration

registry = patronus_registry();

validate_config(versions,registry);

constraint_fcn = registry.constraints(char(versions.constraints));
objective_fcn  = registry.objective(char(versions.objective));



%% Initialize configuration

[params, bounds] = cruise_config(versions.AC_version);

%% Generate model lookups

models = models_config(versions);

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

[c, ceq] = constraint_fcn(X_struct, models, params, bounds);

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
    sol_obj = opti.solve();

    sol_raw     = sol_obj.value(X);
    fval        = sol_obj.value(obj_fun);
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
    [results, checks, aux] = cruise_postp(sol_obj, X, models, params, fval,versions);
    generate_results_page(results, models)
else
     warning('Skipping post-processing: no converged solution.');
end  

% IPOPT Diagnosis 

stats = opti.stats();
inf_pr = stats.iterations.inf_pr;   % primal infeasibility (constraint violation) per iter
inf_du = stats.iterations.inf_du;   % dual infeasibility (KKT stationarity) per iter
obj_hist = stats.iterations.obj;

viol  = computeConstraintViolation(hist.g_hist, n_ineq, n_eq);
slack = computeBoundSlack(hist.X_hist, bounds);

history_matrix = hist.X_hist;
plotIpoptconvergence(history_matrix, bounds);
plotVariableEvolution(hist.X_hist, bounds);
plotConvergenceMetrics(stats, viol, slack);