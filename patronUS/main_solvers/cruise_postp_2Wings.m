function [results, checks, aux] = cruise_postp_2Wings(sol_obj,X ,models, params,fval, versions)
varNames = fieldnames(X);
for i = 1:numel(varNames)
    varName = varNames{i};
    results.opt.(varName) = sol_obj.value(X.(varName));
end

% 2. Calculate the optimized J and phi coordinates
results.J1   = results.opt.V / (results.opt.n1 * params.prop.diameter);
results.phi1 = results.opt.alpha + results.opt.epsilon1;

results.J2   = results.opt.V / (results.opt.n2 * params.prop.diameter);
results.phi2 = results.opt.alpha + results.opt.epsilon2;

% 3. Evaluate models numerically using full()

results.CT1 = full(models.CT_lookup(results.J1,results.phi1));
results.CP1 = full(models.CP_lookup(results.J1,results.phi1));
results.CH1 = full(models.CH_lookup(results.J1,results.phi1));

results.CT2 = full(models.CT_lookup(results.J2,results.phi2));
results.CP2 = full(models.CP_lookup(results.J2,results.phi2));
results.CH2 = full(models.CH_lookup(results.J2,results.phi2));

results.xi_dw =  full (downwash_calc(results.opt.alpha));

results.CL_wing = full(models.CL_wing_lookup(results.opt.alpha));
results.CD_wing = full(models.CD_wing_lookup(results.opt.alpha));
results.CM_wing = full(models.CM_wing_lookup(results.opt.alpha));

results.CL_tail = full(models.CL_tail_lookup(results.opt.alpha - results.xi_dw));
results.CD_tail = full(models.CD_tail_lookup(results.opt.alpha - results.xi_dw));
results.CM_tail = full(models.CM_tail_lookup(results.opt.alpha - results.xi_dw));

results.CL_fuselage = full(models.CL_fuselage_lookup(results.opt.alpha));
results.CD_fuselage = full(models.CD_fuselage_lookup(results.opt.alpha));
results.CM_fuselage = full(models.CM_fuselage_lookup(results.opt.alpha));


% Evaluate result

results.x_range = 1/fval *params.prop.diameter* params.prop.max_rps /(params.prop.P_max_eng) * params.prop.eff * params.E_batt;

%% Calculate derived forces:

q = 0.5 * params.rho * results.opt.V^2; % Dynamic pressure

results.L_wing = results.CL_wing * q * params.wing_area;
results.D_wing = results.CD_wing * q * params.wing_area;
results.M_wing = results.CM_wing * q * params.wing_area*params.geo.c_wing;

results.L_tail = results.CL_tail * q * params.wing_area;
results.D_tail = results.CD_tail * q * params.wing_area;
results.M_tail = results.CM_tail * q * params.wing_area*params.geo.c_wing;

results.L_fuselage = results.CL_fuselage * q * params.fus_area;
results.D_fuselage = results.CD_fuselage * q * params.fus_area;
results.M_fuselage = results.CM_fuselage * q * params.fus_area*params.fus_Lref;


results.W = params.mass * params.g;

% Thrust computation (Assuming T = CT * rho * n^2 * D^4)
results.T1 = params.prop.num_engines* results.CT1 * params.rho * (results.opt.n1)^2 * (params.prop.diameter)^4;
results.H1 = params.prop.num_engines* results.CH1 * params.rho * (results.opt.n1)^2 * (params.prop.diameter)^4;
results.P1 = params.prop.num_engines* results.CP1 * params.rho * (results.opt.n1)^3 * (params.prop.diameter)^5;
aux.propeff1  = results.J1*results.CT1/results.CP1 * 100;

results.T2 = params.prop.num_engines* results.CT2 * params.rho * (results.opt.n2)^2 * (params.prop.diameter)^4;
results.H2 = params.prop.num_engines* results.CH2 * params.rho * (results.opt.n2)^2 * (params.prop.diameter)^4;
results.P2 = params.prop.num_engines* results.CP2 * params.rho * (results.opt.n2)^3 * (params.prop.diameter)^5;
aux.propeff2  = results.J2*results.CT2/results.CP2 * 100;

%% Calculate induced speed

% Nac1
kappa1 = results.T1/(0.5 * params.rho * params.prop.diameter^2*results.opt.V^2);
f_lambda1 = @(lambda_i) lambda_i^4 + 2*sin(results.phi1)*lambda_i^3 + lambda_i^2 - kappa1^2;
lambda_i1 = fzero(f_lambda1 ,0.1);
results.alpha_nac1 = atan(sin(results.phi1)/(lambda_i1 + cos(results.phi1)));
results.V_nac1 = sqrt(results.opt.V^2 + (lambda_i1*results.opt.V)^2 + 2*cos(results.phi1)*lambda_i1*results.opt.V^2);

q_nac1 = 0.5 * params.rho * results.V_nac1^2; % Dynamic pressure

results.CD_nac1 = full(models.CD_fuselage_lookup(results.alpha_nac1));
results.D_nac1 = results.CD_nac1* q_nac1*params.prop.nacelle_area;

results.vi1 = lambda_i1*results.opt.V;
aux.lambda_i1 = lambda_i1;
aux.xi_nac1 = results.phi1-results.alpha_nac1;

% Nac2
kappa2 = results.T2/(0.5 * params.rho * params.prop.diameter^2*results.opt.V^2);
f_lambda2 = @(lambda_i) lambda_i^4 + 2*sin(results.phi2)*lambda_i^3 + lambda_i^2 - kappa2^2;
lambda_i2 = fzero(f_lambda2 ,0.1);
results.alpha_nac2 = atan(sin(results.phi2)/(lambda_i2 + cos(results.phi2)));
results.V_nac2 = sqrt(results.opt.V^2 + (lambda_i2*results.opt.V)^2 + 2*cos(results.phi2)*lambda_i2*results.opt.V^2);

q_nac2 = 0.5 * params.rho * results.V_nac2^2; % Dynamic pressure

results.CD_nac2 = full(models.CD_fuselage_lookup(results.alpha_nac2));
results.D_nac2 = results.CD_nac2* q_nac2*params.prop.nacelle_area;
results.vi2 = lambda_i2*results.opt.V;
aux.lambda_i2 = lambda_i2;
aux.xi_nac2 = results.phi2-results.alpha_nac2;

%% Sum Aero Forces

results.D_Tot = results.D_wing + results.D_fuselage + results.D_tail*cos(results.xi_dw) + ...
    results.D_nac1*cos(results.phi1-results.alpha_nac1) + results.D_nac2*cos(results.phi2-results.alpha_nac2) + results.L_tail*sin(results.xi_dw);
results.L_Tot = results.L_wing + results.L_fuselage + results.L_tail*cos(results.xi_dw) -...
    results.D_nac1*sin(results.phi1-results.alpha_nac1) - results.D_nac2*sin(results.phi2-results.alpha_nac2) -results.D_tail*sin(results.xi_dw);

%% Calculate auxiliary derivative results

aux.AeroEff_wing = results.CL_wing/results.CD_wing;
aux.AeroEff_fuselage = results.CL_fuselage/results.CD_fuselage;
aux.AeroEff_tail = results.CL_tail/results.CD_tail;

aux.AeroEff_Tot = results.L_Tot/results.D_Tot;

%% Sum of aero Forces per component in wind axes ( for moments eq)

results.L_tailwing_wind = results.L_tail*cos(results.xi_dw) - results.D_tail*sin(results.xi_dw);
results.D_tailwing_wind = results.D_tail*cos(results.xi_dw) + results.L_tail*sin(results.xi_dw);

%% Definition of Terms for the moments equation ( to ease up the constraint shape)

results.MA_wing = q*params.wing_area*params.geo.c_wing*results.CM_wing; % Aero moment
results.MA_fus  = q*params.fus_area*params.fus_Lref*results.CM_fuselage; % Aero moment
results.MA_tail = q*params.wing_area*params.geo.c_wing*results.CM_tail; % Aero moment

results.MA = results.MA_wing + results.MA_fus + results.MA_tail;

% Inertia_nac1 = params.I_yy1 *epsilon1dotdot;
% Inertia_nac2 = params.I_yy2 *epsilon2dotdot;

results.M_Fus = results.L_fuselage*cos(results.opt.alpha)*params.geo.xfus - results.L_fuselage*sin(results.opt.alpha)*params.geo.zfus...
    +results.D_fuselage*sin(results.opt.alpha)*params.geo.xfus - results.D_fuselage*cos(results.opt.alpha)*params.geo.zfus; % Fuselage contribution

results.M_Wing = results.L_wing* cos(results.opt.alpha)*params.geo.xw -results.L_wing*sin(results.opt.alpha)*params.geo.zw...
    +results.D_wing*sin(results.opt.alpha)*params.geo.xw + results.D_wing*cos(results.opt.alpha)*params.geo.zw; % Wing contribution

results.M_Tail = - results.L_tailwing_wind*cos(results.opt.alpha)*params.geo.xtw - results.L_tailwing_wind*sin(results.opt.alpha)*params.geo.ztw...
    -results.D_tailwing_wind*sin(results.opt.alpha)*params.geo.xtw + results.D_tailwing_wind*cos(results.opt.alpha)*params.geo.ztw; % Tail wing contribution

results.M_Eng1 = results.T1*sin(results.opt.epsilon1)*params.geo.xtw - results.T1*cos(results.opt.epsilon1)*params.geo.ztw ...
    + results.H1*cos(results.opt.epsilon1)*params.geo.xtw + results.H1*sin(results.opt.epsilon1)*params.geo.ztw;

results.M_Eng2 = -results.T2*sin(results.opt.epsilon2)*params.geo.xtw - results.T2*cos(results.opt.epsilon2)*params.geo.ztw ...
    - results.H2*cos(results.opt.epsilon2)*params.geo.xtw + results.H2*sin(results.opt.epsilon2)*params.geo.ztw;

results.M_Nac1 = results.D_nac1*sin(results.opt.epsilon1-results.alpha_nac1)*params.geo.xw + ...
    results.D_nac1*cos(results.opt.epsilon1-results.alpha_nac1)*params.geo.zw;

results.M_Nac2 = -results.D_nac2*sin(results.opt.epsilon2-results.alpha_nac2)*params.geo.xtw + ...
    results.D_nac2*cos(results.opt.epsilon2-results.alpha_nac2)*params.geo.ztw;


% Evaluate Check in force balance

checks.horcheck = results.T1*cos(results.phi1) + results.T2*cos(results.phi2) - results.H1*sin(results.phi1)- results.H2*sin(results.phi2)-results.D_Tot;
checks.vercheck = results.T1*sin(results.phi1) + results.T2*sin(results.phi2) + results.H1*cos(results.phi1)+ results.H2*cos(results.phi2)+results.L_Tot...
    -params.mass*params.g;
checks.momcheck = results.MA + results.M_Wing + results.M_Fus + results.M_Tail+  results.M_Nac1 + results.M_Nac2 + results.M_Eng1 + results.M_Eng2;