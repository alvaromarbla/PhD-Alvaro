function F = vtolto_forces_3GDL(Z, models, params)
% Forces and pitching moment of the 2Wings aircraft at one trajectory node.
% Same aero, propulsive and nacelle models as cruise_constraints_3GDL, but the
% resultant is written in inertial axes (x forward, h up) so it stays defined
% in hover, where gamma and the wind axes are not.

vx       = Z.vx;
vh       = Z.vh;
theta    = Z.theta;
epsilon1 = Z.epsilon1;
epsilon2 = Z.epsilon2;
n1       = Z.n1;
n2       = Z.n2;

%% Relative wind (regularised with V_eps so everything is smooth at hover)

V  = sqrt(vx^2 + vh^2 + params.V_eps^2);
cg = vx/V;  % cos(gamma)
sg = vh/V;  % sin(gamma)

% Body axes velocities. alpha = theta - gamma, -90 deg in a vertical climb.
% Valid for u_b >= 0 (no backwards flight)
u_b = vx*cos(theta) + vh*sin(theta);
w_b = vx*sin(theta) - vh*cos(theta);
alpha = atan2(w_b, sqrt(u_b^2 + params.V_eps^2));
gamma = theta - alpha;

xi_dw = downwash_calc(alpha);

%% Propulsive coefficients

J1   = V/(n1*params.prop.diameter);
phi1 = alpha + epsilon1;
CT1 = models.CT_lookup(J1,phi1);
CP1 = models.CP_lookup(J1,phi1);
CH1 = models.CH_lookup(J1,phi1);

J2   = V/(n2*params.prop.diameter);
phi2 = alpha + epsilon2;
CT2 = models.CT_lookup(J2,phi2);
CP2 = models.CP_lookup(J2,phi2);
CH2 = models.CH_lookup(J2,phi2);

%% Aero coefficients

CL_wing = models.CL_wing_lookup(alpha);
CD_wing = models.CD_wing_lookup(alpha);
CM_wing = models.CM_wing_lookup(alpha);

CL_fus = models.CL_fuselage_lookup(alpha);
CD_fus = models.CD_fuselage_lookup(alpha);
CM_fus = models.CM_fuselage_lookup(alpha);

CL_tailwing = models.CL_tail_lookup(alpha-xi_dw);
CD_tailwing = models.CD_tail_lookup(alpha-xi_dw);
CM_tailwing = models.CM_tail_lookup(alpha-xi_dw);

%% Aero forces (wind axes)

q = 0.5 * params.rho * (V^2); % Dynamic pressure

L_wing = CL_wing*q* params.wing_area;
D_wing = CD_wing*q* params.wing_area;

L_fus = CL_fus*q* params.fus_area;
D_fus = CD_fus*q* params.fus_area;

L_tailwing = CL_tailwing*q* params.tailwing_area;
D_tailwing = CD_tailwing*q* params.tailwing_area;

L_tailwing_wind = L_tailwing*cos(xi_dw) - D_tailwing*sin(xi_dw);
D_tailwing_wind = D_tailwing*cos(xi_dw) + L_tailwing*sin(xi_dw);

L_aero = L_wing + L_fus + L_tailwing_wind;
D_aero = D_wing + D_fus + D_tailwing_wind;

%% Prop forces

T1 =  params.prop.num_engines * CT1 * params.rho * (n1^2) * (params.prop.diameter^4);
H1 =  params.prop.num_engines * CH1 * params.rho * (n1^2) * (params.prop.diameter^4);
P1 =  params.prop.num_engines * CP1 * params.rho * (n1^3) * (params.prop.diameter^5);

T2 =  params.prop.num_engines * CT2 * params.rho * (n2^2) * (params.prop.diameter^4);
H2 =  params.prop.num_engines * CH2 * params.rho * (n2^2) * (params.prop.diameter^4);
P2 =  params.prop.num_engines * CP2 * params.rho * (n2^3) * (params.prop.diameter^5);

%% Nacelle drag

[D_nac1, alpha_nac1, vi1] = nacelle_drag(T1, phi1, V, models, params);
[D_nac2, alpha_nac2, vi2] = nacelle_drag(T2, phi2, V, models, params);

%% Resultant in inertial axes
% Aero: wind axes rotated by gamma. Thrust/H: angle theta+epsilon from the
% horizon. Nacelle drag: angle theta+epsilon-alpha_nac (as in cruise for gamma = 0)

Fx = - D_aero*cg - L_aero*sg ...
     + T1*cos(theta+epsilon1) - H1*sin(theta+epsilon1) ...
     + T2*cos(theta+epsilon2) - H2*sin(theta+epsilon2) ...
     - D_nac1*cos(theta+epsilon1-alpha_nac1) - D_nac2*cos(theta+epsilon2-alpha_nac2);

Fh = - D_aero*sg + L_aero*cg ...
     + T1*sin(theta+epsilon1) + H1*cos(theta+epsilon1) ...
     + T2*sin(theta+epsilon2) + H2*cos(theta+epsilon2) ...
     - D_nac1*sin(theta+epsilon1-alpha_nac1) - D_nac2*sin(theta+epsilon2-alpha_nac2) ...
     - params.mass*params.g;

%% Pitching moment (same terms as cruise_constraints_3GDL)

MA_wing = q*params.wing_area*params.geo.c_wing*CM_wing;
MA_fus  = q*params.fus_area*params.fus_Lref*CM_fus;
MA_tail = q*params.wing_area*params.geo.c_wing*CM_tailwing;

M_Fus = L_fus*cos(alpha)*params.geo.xfus - L_fus*sin(alpha)*params.geo.zfus...
    +D_fus*sin(alpha)*params.geo.xfus  - D_fus*cos(alpha)*params.geo.zfus;

M_Wing = L_wing* cos(alpha)*params.geo.xw -L_wing*sin(alpha)*params.geo.zw...
    +D_wing*sin(alpha)*params.geo.xw + D_wing*cos(alpha)*params.geo.zw;

M_Tail = - L_tailwing_wind*cos(alpha)*params.geo.xtw - L_tailwing_wind*sin(alpha)*params.geo.ztw...
    -D_tailwing_wind*sin(alpha)*params.geo.xtw + D_tailwing_wind*cos(alpha)*params.geo.ztw;

M_Eng1 = T1*sin(epsilon1)*params.geo.xtw - T1*cos(epsilon1)*params.geo.ztw ...
    + H1*cos(epsilon1)*params.geo.xtw + H1*sin(epsilon1)*params.geo.ztw;

M_Eng2 = -T2*sin(epsilon2)*params.geo.xtw - T2*cos(epsilon2)*params.geo.ztw ...
    - H2*cos(epsilon2)*params.geo.xtw + H2*sin(epsilon2)*params.geo.ztw;

M_Nac1 = D_nac1*sin(epsilon1-alpha_nac1)*params.geo.xw + D_nac1*cos(epsilon1-alpha_nac1)*params.geo.zw;

M_Nac2 = -D_nac2*sin(epsilon2-alpha_nac2)*params.geo.xtw + D_nac2*cos(epsilon2-alpha_nac2)*params.geo.ztw;

M = MA_wing + MA_fus + MA_tail + M_Wing + M_Fus + M_Tail + M_Nac1 + M_Nac2 + M_Eng1 + M_Eng2;

%% Output

F.Fx = Fx;      F.Fh = Fh;      F.M = M;
F.V = V;        F.alpha = alpha; F.gamma = gamma; F.xi_dw = xi_dw;
F.phi1 = phi1;  F.phi2 = phi2;
F.J1 = J1;      F.J2 = J2;
F.CT1 = CT1;    F.CT2 = CT2;
F.CP1 = CP1;    F.CP2 = CP2;
F.CH1 = CH1;    F.CH2 = CH2;
F.T1 = T1;      F.T2 = T2;
F.H1 = H1;      F.H2 = H2;
F.P1 = P1;      F.P2 = P2;
F.CL_wing = CL_wing; F.CD_wing = CD_wing; F.CM_wing = CM_wing;
F.CL_tail = CL_tailwing; F.CD_tail = CD_tailwing; F.CM_tail = CM_tailwing;
F.CL_fuselage = CL_fus; F.CD_fuselage = CD_fus; F.CM_fuselage = CM_fus;
F.L_aero = L_aero; F.D_aero = D_aero;
F.D_nac1 = D_nac1; F.D_nac2 = D_nac2;
F.alpha_nac1 = alpha_nac1; F.alpha_nac2 = alpha_nac2;
F.vi1 = vi1;    F.vi2 = vi2;
F.MA = MA_wing + MA_fus + MA_tail;
F.M_Wing = M_Wing; F.M_Fus = M_Fus; F.M_Tail = M_Tail;
F.M_Eng1 = M_Eng1; F.M_Eng2 = M_Eng2;
F.M_Nac1 = M_Nac1; F.M_Nac2 = M_Nac2;

end

function [D_nac, alpha_nac, vi] = nacelle_drag(T, phi, V, models, params)
% Cruise nacelle model multiplied through by V^4 (lambda_i = vi/V), so it
% stays finite at hover:  vi^4 + 2*sin(phi)*V*vi^3 + V^2*vi^2 - K^2 = 0
% Newton starts from the hover solution vi = sqrt(K)

K  = T/(0.5 * params.rho * params.prop.diameter^2);
vi = sqrt(sqrt(K^2 + 1e-6));

for i = 1:8
    f  = vi^4 + 2*sin(phi)*V*vi^3 + V^2*vi^2 - K^2;
    df = 4*vi^3 + 6*sin(phi)*V*vi^2 + 2*V^2*vi;
    vi = vi - f/df;
end

alpha_nac = atan2(V*sin(phi), vi + V*cos(phi));
CD_nac    = models.CD_nac_lookup(alpha_nac);

V_nac2 = V^2 + vi^2 + 2*cos(phi)*vi*V;
q_nac  = 0.5 * params.rho * V_nac2;

D_nac = CD_nac*q_nac* params.prop.nacelle_area*params.prop.num_engines; % Area accounts for ONE engine

end
