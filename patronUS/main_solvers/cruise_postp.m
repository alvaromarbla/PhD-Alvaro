function [results, checks, aux] = cruise_postp(sol_obj,X ,models, params,fval, versions)

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

results.x_range_opt = 1/fval *params.prop.diameter* params.prop.max_rps /(params.prop.P_max_eng) * params.prop.eff * params.E_batt;

%% Calculate auxiliary derivative results

aux.AeroEff = results.CL_opt/results.CD_opt;
%% Calculate derived forces:

q = 0.5 * params.rho * results.V_opt^2; % Dynamic pressure
results.L_opt = results.CL_opt * q * params.wing_area;
results.D_opt = results.CD_opt * q * params.wing_area;
results.W_opt = params.mass * params.g;

% Thrust computation (Assuming T = CT * rho * n^2 * D^4)
results.T_opt = params.prop.num_engines* results.CT_opt * params.rho * (results.n_opt)^2 * (params.prop.diameter)^4;
results.P_opt = params.prop.num_engines* results.CP_opt * params.rho * (results.n_opt)^3 * (params.prop.diameter)^5;
aux.propeff  = results.J_opt*results.CT_opt/results.CP_opt * 100;

%% Calculate induced speed

if versions.constraints == "1Wing_Nacelle"
    kappa = results.T_opt/(0.5 * params.rho * params.prop.diameter^2*results.V_opt^2);
    f_lambda = @(lambda_i) lambda_i^4 + 2*sin(results.phi_opt)*lambda_i^3 + lambda_i^2 - kappa^2;
    lambda_i = fzero(f_lambda ,0.1);
    results.alpha_nac = atan(sin(results.phi_opt)/(lambda_i + cos(results.phi_opt)));
    results.V_nac = sqrt(results.V_opt^2 + (lambda_i*results.V_opt)^2 + 2*cos(results.phi_opt)*lambda_i*results.V_opt^2);

    q_nac = 0.5 * params.rho * results.V_nac^2; % Dynamic pressure

    results.CD_nac_opt = full(models.CD_fuselage_lookup(results.alpha_nac));
    results.D_nac_opt = results.CD_nac_opt* q_nac*params.nacelle_area;
    results.vi = lambda_i*results.V_opt;
    aux.lambda_i = lambda_i;
    aux.xi_nac = results.phi_opt-results.alpha_nac;
end

% Evaluate Check in force balance

checks.horcheck = results.T_opt*cos(results.phi_opt) - results.D_opt - results.D_nac_opt*cos(results.phi_opt - results.alpha_nac);
checks.vercheck = results.L_opt + results.T_opt*sin(results.phi_opt)- results.D_nac_opt*sin(results.phi_opt - results.alpha_nac) - (params.mass * params.g);