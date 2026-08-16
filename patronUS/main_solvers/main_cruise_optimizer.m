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
versions.CL_version   = 'dummy'; % CL model to use ("analytical, dummy")
versions.CD_version   = 'dummy'; % CD model to use ("analytical, dummy")
versions.CM_version   = 'dummy'; % CM model to use ("analytical, dummy")


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

[c, ceq, diag] = constraint_fcn(X_struct, models, params, bounds); % Can have diagnosis ( or not) 


%%%%%%%%%%%%%%%%%% EVALUATE INITIAL ( TRIM ) POINT %%%%%%%%%%%%%%%%%%

if exist("diag","var")

diag_names = fieldnames(diag);
    diag_exprs = {};
    for k = 1:numel(diag_names)
        diag_exprs{end+1} = diag.(diag_names{k});
    end

    % Build a CasADi Function: X (numeric vector) -> all diag fields
    trim_eval = Function('trim_eval', {X}, diag_exprs, {'X'}, diag_names);

    % Evaluate at the initial guess vars0 (plain doubles)
    diag_vals_cell = cell(1, numel(diag_names));
    [diag_vals_cell{:}] = trim_eval(vars0);

    % Repack into a plain MATLAB struct of doubles
    diag0 = struct();
    for k = 1:numel(diag_names)
        diag0.(diag_names{k}) = full(diag_vals_cell{k});
    end

    fprintf('alpha = %.3f deg | phi1 = %.3f deg | phi2 = %.3f deg\n', ...
        rad2deg(diag0.alpha), rad2deg(diag0.phi1), rad2deg(diag0.phi2));
    fprintf('n1 = %.3f rps | n2 = %.3f rps\n', ...
        diag0.n1, diag0.n2);


    fprintf('\n--- FORCE BREAKDOWN [N] ---\n');
    fprintf('T1 = %10.2f   T2 = %10.2f   H1 = %10.2f   H2 = %10.2f\n', diag0.T1, diag0.T2, diag0.H1, diag0.H2);
    fprintf('D_wing = %8.2f  D_fus = %8.2f  D_tail = %8.2f  D_nac1 = %8.2f  D_nac2 = %8.2f  D_Tot = %8.2f\n', ...
        diag0.D_wing, diag0.D_fus, diag0.D_tailwing, diag0.D_nac1, diag0.D_nac2, diag0.D_Tot);
    fprintf('L_wing = %8.2f  L_fus = %8.2f  L_tail = %8.2f  L_Tot  = %8.2f  Weight = %8.2f\n', ...
        diag0.L_wing, diag0.L_fus, diag0.L_tailwing, diag0.L_Tot, diag0.weight);

    fprintf('\n--- MOMENT BREAKDOWN [N*m] ---\n');
    fprintf('MA_wing = %8.2f  MA_fus = %8.2f  MA_tail = %8.2f\n', diag0.MA_wing, diag0.MA_fus, diag0.MA_tail);
    fprintf('M_Wing  = %8.2f  M_Fus  = %8.2f  M_Tail  = %8.2f\n', diag0.M_Wing, diag0.M_Fus, diag0.M_Tail);
    fprintf('M_Eng1  = %8.2f  M_Eng2 = %8.2f  M_Nac1  = %8.2f  M_Nac2 = %8.2f\n', ...
        diag0.M_Eng1, diag0.M_Eng2, diag0.M_Nac1, diag0.M_Nac2);

    fprintf('\n--- CEQ RESIDUALS ---\n');
    fprintf('Long force (should be 0): %10.3f  N\n', diag0.ceq(1));
    fprintf('Trans force (should be 0): %10.3f  N\n', diag0.ceq(2));
    fprintf('Moment (should be 0):      %10.3f  N*m\n', diag0.ceq(3));
end


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