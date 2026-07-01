function [results, aux] = cruise_postp(sol_obj,X ,models, params,fval)

%% Postprocess
% 1. Extract the optimized numerical values
results.V_opt       = sol_obj.value(X(1));
results.gamma_opt   = sol_obj.value(X(2));
results.alpha_opt   = sol_obj.value(X(3));
results.epsilon_opt = sol_obj.value(X(4));
results.n_opt       = sol_obj.value(X(5));

% 2. Calculate the optimized J and phi coordinates
results.J_opt   = results.V_opt / (results.n_opt * params.prop.diameter);
results.phi_opt = results.alpha_opt + results.epsilon_opt;

% 3. Evaluate models numerically using full()
results.CT_opt = full(models.CT_lookup(results.J_opt,results.phi_opt));
results.CP_opt = full(models.CP_lookup(results.J_opt,results.phi_opt));
results.CL_opt = full(models.CL_lookup(results.alpha_opt));
results.CD_opt = full(models.CD_lookup(results.alpha_opt));

% Evaluate result

results.x_range_opt = 1/fval;

%% Calculate auxiliary derivative results

aux.AeroEff = results.CL_opt/results.CD_opt;

%% Calculate derived forces:

q = 0.5 * params.rho * results.V_opt^2; % Dynamic pressure
results.L_opt = results.CL_opt * q * params.wing_area;
results.D_opt = results.CD_opt * q * params.wing_area;
results.W_opt = params.mass * params.g;

% Thrust computation (Assuming T = CT * rho * n^2 * D^4)
results.T_opt = results.CT_opt * params.rho * (results.n_opt)^2 * (params.prop.diameter)^4;

%% Call forces plot calculator

plotforceangles(results);