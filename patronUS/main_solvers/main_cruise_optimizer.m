%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PATRONUS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Performance Analysis, Trajectory Routing
% and Optimized Navigation for Unmanned Systems

%% Main Cruise Segment Optimizer
% By: Alvaro Martinez Blanco (2026)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User version input

versions.AC_version   = "2WING"; % Aicraft Version
versions.CL_version   = 'analytical'; % CL model to use
versions.CD_version   = 'analytical'; % CD model to use
versions.CT_version   = 'analytical'; % CT model to use
versions.CP_version   = 'analytical'; % CP model to use
versions.Batt_version = ' '; % Battery model to use

% Versions for constraints
% To choose among: 1Wing, 1Wing_Nacelle, 2Wings
versions.constraints = "2Wings";

%% Import Casadi and optimizer

import casadi.*
opti = Opti();

%% Initialize configuration

[params, bounds] = cruise_config(versions.AC_version);

%% Generate model lookups

models = models_config(versions);

%% Define Decision Variables

X = opti.variable(8, 1); % [V; gamma; alpha; epsilon1 ; epsilon2; n1; n2;delta_e]

% Apply bounds to states

opti.subject_to(bounds.V.min        <= X(1) <= bounds.V.max);
opti.subject_to(bounds.gamma.min    <= X(2) <= bounds.gamma.max);
opti.subject_to(bounds.alpha.min    <= X(3) <= bounds.alpha.max);
opti.subject_to(bounds.epsilon1.min  <= X(4) <= bounds.epsilon1.max);
opti.subject_to(bounds.n1.min        <= X(5) <= bounds.n1.max);
opti.subject_to(bounds.epsilon2.min  <= X(6) <= bounds.epsilon2.max);
opti.subject_to(bounds.n2.min        <= X(7) <= bounds.n2.max);
opti.subject_to(bounds.delta_e.min   <= X(8) <= bounds.delta_e.max);

%% Apply initial interation for solver

vars0 = [24.225109; deg2rad(0); 0.046444; deg2rad(60); 25.214554; deg2rad(60); 25.214554 ; deg2rad(0)];
%vars0 = [30; deg2rad(3); 0.1; deg2rad(4); 40];
opti.set_initial(X, vars0)

%% Call Objective Function

obj_fun = cruise_objective(X, models.CP_lookup, params);
opti.minimize(obj_fun);

%% Call Constraints
switch versions.constraints
    case "1Wing"
        [c, ceq] = cruise_constraints_2GDL(X, models, params,bounds);
    case "1Wing_Nacelle"
        [c, ceq] = cruise_constraints_nacelle(X, models, params,bounds);
    case "2Wings"
        [c, ceq] = cruise_constraints_3GDL(X, models, params,bounds);
    otherwise
        error("Unknown constraint set selected.")
end

g_all = [c; ceq];   % keep this handle
n_ineq = numel(c);
n_eq   = numel(ceq);

opti.subject_to(c <= 0);
opti.subject_to(ceq == 0);

hist = IterHistory();
opti.callback(@(i) hist.record(opti.debug.value(X), opti.debug.value(g_all)));
%% Configure solver (IPOPT)

opts.ipopt.max_iter = 1000;
opts.ipopt.tol = 1e-9; % ConstraintTolerance
opti.solver('ipopt', opts);

%% Run

try
    sol_obj = opti.solve();

    % 4. Extract numerical results (Equivalent to your [sol, fval])
    sol  = sol_obj.value(X);     % Optimized decision variables vector
    fval = sol_obj.value(obj_fun); % Optimized objective function value

catch e
    fprintf('Optimization failed or hit limits. Extracting last debug values.\n');
    sol  = opti.debug.value(X);
    fval = opti.debug.value(obj_fun);
end

[results, aux] = cruise_postp(sol_obj, X , models, params, fval);

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