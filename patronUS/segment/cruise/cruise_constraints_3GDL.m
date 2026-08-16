function [c, ceq] = cruise_constraints_3GDL(X_struct, models, params, bounds)
% Steady State constraints

V       = X_struct.V;
gamma   = X_struct.gamma; % Steady cruise implies gamma = 0, which is handled via bounds or optimization rules
theta   = X_struct.theta;
epsilon1 = X_struct.epsilon1;
epsilon2 = X_struct.epsilon2;
n1      = X_struct.n1;
n2      = X_struct.n2;
% tau1     = x(8);
% tau2     = x(9);
delta_e = X_struct.deltae;

alpha = theta-gamma;
xi_dw = downwash_calc(alpha); %%%%%%%

%% Call B-spline models on propulsive forces:
J1 = V./(n1*params.prop.diameter);
phi1 = alpha + epsilon1;

CT1 = models.CT_lookup(J1,phi1);
CP1 = models.CP_lookup(J1,phi1);
CH1 = models.CH_lookup(J1,phi1);

J2 = V./(n2*params.prop.diameter);
phi2 = alpha + epsilon2;

CT2 = models.CT_lookup(J2,phi2);
CP2 = models.CP_lookup(J2,phi2);
CH2 = models.CH_lookup(J2,phi2);

%% Call B-spline models on aero forces:

CL_wing = models.CL_wing_lookup(alpha);
CD_wing = models.CD_wing_lookup(alpha);
CM_wing = models.CM_wing_lookup(alpha);

CL_fus = models.CL_fuselage_lookup(alpha);
CD_fus = models.CD_fuselage_lookup(alpha);
CM_fus = models.CM_fuselage_lookup(alpha);

CL_tailwing = models.CL_tail_lookup(alpha-xi_dw); %, delta_e)
CD_tailwing = models.CD_tail_lookup(alpha-xi_dw); %, delta_e)
CM_tailwing = models.CM_tail_lookup(alpha); %, delta_e)

%% Force calculations

q = 0.5 * params.rho * (V^2); % Dynamic pressure

L_wing = CL_wing*q* params.wing_area;
D_wing = CD_wing*q* params.wing_area;

L_fus = CL_fus*q* params.fus_area;
D_fus = CD_fus*q* params.fus_area;

L_tailwing = CL_tailwing*q* params.tailwing_area;
D_tailwing = CD_tailwing*q* params.tailwing_area;

%% Prop forces calculation
T1 =  params.prop.num_engines * CT1 * params.rho * (n1^2) * (params.prop.diameter^4);
H1 =  params.prop.num_engines * CH1 * params.rho * (n1^2) * (params.prop.diameter^4);
P1 =  params.prop.num_engines * CP1 * params.rho * (n1^3) * (params.prop.diameter^5);

T2 =  params.prop.num_engines * CT2 * params.rho * (n2^2) * (params.prop.diameter^4);
H2 =  params.prop.num_engines * CH2 * params.rho * (n2^2) * (params.prop.diameter^4);
P2 =  params.prop.num_engines * CP2 * params.rho * (n2^3) * (params.prop.diameter^5);

%% Nacelle drag Calculation

% NACELLE 1
kappa1 = T1/(0.5 * params.rho * params.prop.diameter^2*V^2);

lambda_i1 = 0.1;

for i= 1:5
    f = lambda_i1^4 + 2*sin(phi1)*lambda_i1^3 + lambda_i1^2 - kappa1^2;

    df = 4*lambda_i1^3 + 6*sin(phi1)*lambda_i1^2 + 2*lambda_i1;

    lambda_i1 = lambda_i1 -f/df;
end

% Filter for the real, positive physical root corresponding to induced velocity
alpha_nac1 = atan(sin(phi1)/(lambda_i1 + cos(phi1)));
CD_nac1 = models.CD_nac_lookup(alpha_nac1);


V1_nac2 = V^2 + (lambda_i1*V)^2 + 2*cos(phi1)*lambda_i1*V^2;
q_nac1 = 0.5 * params.rho * (V1_nac2) ;

D_nac1 = CD_nac1*q_nac1* params.prop.nacelle_area*params.prop.num_engines; % Area accounts for ONE engine

% NACELLE 2
kappa2 = T2/(0.5 * params.rho * params.prop.diameter^2*V^2);

lambda_i2 = 0.1;

for i= 1:5
    f = lambda_i2^4 + 2*sin(phi2)*lambda_i2^3 + lambda_i2^2 - kappa2^2;

    df = 4*lambda_i2^3 + 6*sin(phi2)*lambda_i2^2 + 2*lambda_i2;

    lambda_i2 = lambda_i2 -f/df;
end

% Filter for the real, positive physical root corresponding to induced velocity
alpha_nac2 = atan(sin(phi2)/(lambda_i2 + cos(phi2)));
CD_nac2 = models.CD_nac_lookup(alpha_nac2);


V2_nac2 = V^2 + (lambda_i2*V)^2 + 2*cos(phi2)*lambda_i2*V^2;
q_nac2 = 0.5 * params.rho * (V2_nac2) ;

D_nac2 = CD_nac2*q_nac2* params.prop.nacelle_area*params.prop.num_engines; % Area accounts for ONE engine

%% Sum aero Forces

D_Tot = D_wing + D_fus + D_tailwing*cos(xi_dw) + D_nac1*cos(phi1-alpha_nac1) + D_nac2*cos(phi2-alpha_nac2) + L_tailwing*sin(xi_dw);
L_Tot = L_wing + L_fus + L_tailwing*cos(xi_dw) - D_nac1*sin(phi1-alpha_nac1) - D_nac2*sin(phi2-alpha_nac2) - D_tailwing*sin(xi_dw);

%% Sum of aero Forces per component in wind axes ( for moments eq)

L_tailwing_wind = L_tailwing*cos(xi_dw) - D_tailwing*sin(xi_dw);
D_tailwing_wind = D_tailwing*cos(xi_dw) + L_tailwing*sin(xi_dw);

%% Definition of Terms for the moments equation ( to ease up the constraint shape)

MA_wing = q*params.wing_area*params.geo.c_wing*CM_wing; % Aero moment
MA_fus  = q*params.fus_area*params.fus_Lref*CM_fus; % Aero moment
MA_tail = q*params.wing_area*params.geo.c_wing*CM_tailwing; % Aero moment

MA = MA_wing + MA_fus + MA_tail;

% Inertia_nac1 = params.I_yy1 *epsilon1dotdot;
% Inertia_nac2 = params.I_yy2 *epsilon2dotdot;

M_Fus = L_fus*cos(alpha)*params.geo.xfus - L_fus*sin(alpha)*params.geo.zfus...
    +D_fus*sin(alpha)*params.geo.xfus  - D_fus*cos(alpha)*params.geo.zfus; % Fuselage contribution

M_Wing = L_wing* cos(alpha)*params.geo.xw -L_wing*sin(alpha)*params.geo.zw...
    +D_wing*sin(alpha)*params.geo.xw + D_wing*cos(alpha)*params.geo.zw; % Wing contribution

M_Tail = - L_tailwing_wind*cos(alpha)*params.geo.xtw - L_tailwing_wind*sin(alpha)*params.geo.ztw...
    -D_tailwing_wind*sin(alpha)*params.geo.xtw + D_tailwing_wind*cos(alpha)*params.geo.ztw; % Tail wing contribution

M_Eng1 = T1*sin(epsilon1)*params.geo.xtw - T1*cos(epsilon1)*params.geo.ztw ...
    + H1*cos(epsilon1)*params.geo.xtw + H1*sin(epsilon1)*params.geo.ztw;

M_Eng2 = -T2*sin(epsilon2)*params.geo.xtw - T2*cos(epsilon2)*params.geo.ztw ...
    - H2*cos(epsilon2)*params.geo.xtw + H2*sin(epsilon2)*params.geo.ztw;

M_Nac1 = D_nac1*sin(epsilon1-alpha_nac1)*params.geo.xw + D_nac1*cos(epsilon1-alpha_nac1)*params.geo.zw;

M_Nac2 = -D_nac2*sin(epsilon2-alpha_nac2)*params.geo.xtw + D_nac2*cos(epsilon2-alpha_nac2)*params.geo.ztw;

%% Constraint satisfaction

% Equalities

ceq = [  T1*cos(phi1) + T2*cos(phi2) - H1*sin(phi1)- H2*sin(phi2)-D_Tot-params.mass*params.g*sin(gamma); % Long forces (m*vdot = ...)
         T1*sin(phi1) + T2*sin(phi2) + H1*cos(phi1)+ H2*cos(phi2)+L_Tot-params.mass*params.g*cos(gamma); % Trans forces (m*V*gammadot = ...)
         MA + M_Wing + M_Fus + M_Tail+  M_Nac1 + M_Nac2 + M_Eng1 + M_Eng2; % Moments (qdot* Iyy = ...)
         gamma % force gamma = 0
    ];

%tau1 = Inertia_nac1*epsilon1dotdot; % Reaction torque for Nac 1
%tau2 = Inertia_nac2*epsilon2dotdot; % Reaction torque for Nac 2

% Inequalities

c = [-V/(n1*params.prop.diameter)  ; ...
    -V/(n2*params.prop.diameter);...
    -CT1;...
    -CT2;...
    -CP1;  ...
    -CP2;  ...
    T1 - params.prop.T_max_eng*params.prop.num_engines;
    T2 - params.prop.T_max_eng*params.prop.num_engines; ...
    P1 - params.prop.P_max_eng*params.prop.num_engines;...
    P2 - params.prop.P_max_eng*params.prop.num_engines;...
    J1-1.2;...
    J2-1.2];

end